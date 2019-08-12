"""
This module contains some tests for Badread. To run them, execute `python3 -m unittest` from the
root Badread directory.

Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Badread

This file is part of Badread. Badread is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Badread is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Badread.
If not, see <http://www.gnu.org/licenses/>.
"""

import collections
import os
import statistics
import unittest
import zlib

import badread.simulate
import badread.fragment_lengths
import badread.misc


class TestLinearFragments(unittest.TestCase):

    def setUp(self):
        self.null = open(os.devnull, 'w')
        self.lengths = badread.fragment_lengths.FragmentLengths(1000, 0, self.null)
        self.ref_seqs = {'r': badread.misc.get_random_sequence(10000)}
        self.ref_depths = {'r': 1.0}
        self.ref_circular = {'r': False}
        self.rev_comp_ref_seqs = {name: badread.misc.reverse_complement(seq)
                                  for name, seq in self.ref_seqs.items()}
        self.ref_contigs, self.ref_contig_weights = \
            badread.simulate.get_ref_contig_weights(self.ref_seqs, self.ref_depths)
        self.trials = 100

    def tearDown(self):
        self.null.close()

    def test_plain(self):
        # No adapters, glitches, chimeras.
        start_adapt_rate, start_adapt_amount = 0.0, 0.0
        end_adapt_rate, end_adapt_amount = 0.0, 0.0
        Args = collections.namedtuple('Args', ['start_adapter_seq', 'end_adapter_seq',
                                               'junk_reads', 'random_reads', 'chimeras',
                                               'glitch_rate', 'glitch_size', 'glitch_skip'])
        args = Args(start_adapter_seq='', end_adapter_seq='',
                    junk_reads=0, random_reads=0, chimeras=0,
                    glitch_rate=0, glitch_size=0, glitch_skip=0)
        forward_count, reverse_count = 0, 0
        lengths = []
        forward_ref, reverse_ref = self.ref_seqs['r'], self.rev_comp_ref_seqs['r']
        for _ in range(self.trials):
            fragment, info = \
                badread.simulate.build_fragment(self.lengths, self.ref_seqs, self.rev_comp_ref_seqs,
                                                self.ref_contigs, self.ref_contig_weights,
                                                self.ref_circular, args, start_adapt_rate,
                                                start_adapt_amount, end_adapt_rate,
                                                end_adapt_amount)
            if fragment in forward_ref:
                forward_count += 1
            if fragment in reverse_ref:
                reverse_count += 1

            # With a linear reference, fragments will either be 1000 bp (the target length) or
            # less (if they land too close to the sequence end).
            self.assertLessEqual(len(fragment), 1000)
            lengths.append(len(fragment))

            # Not full-length fragments should be at the start or end.
            if len(fragment) < 1000:
                self.assertTrue(forward_ref.startswith(fragment) or
                                forward_ref.endswith(fragment) or
                                reverse_ref.startswith(fragment) or
                                reverse_ref.endswith(fragment))

        # Should have some from each strand.
        self.assertGreater(forward_count, 10)
        self.assertGreater(reverse_count, 10)
        self.assertGreaterEqual(forward_count + reverse_count, self.trials)

        # While the average will be less than 1000 bp, 1000 should be the most common length.
        most_common_length = max(set(lengths), key=lengths.count)
        self.assertEqual(most_common_length, 1000)

    def test_full_adapters(self):
        # Start and end adapters reliably put on entirely.
        start_adapt_rate, start_adapt_amount = 1.0, 1.0
        end_adapt_rate, end_adapt_amount = 1.0, 1.0
        Args = collections.namedtuple('Args', ['start_adapter_seq', 'end_adapter_seq',
                                               'junk_reads', 'random_reads', 'chimeras',
                                               'glitch_rate', 'glitch_size', 'glitch_skip'])
        args = Args(start_adapter_seq='TGTATAATACGACGCCGAGC',
                    end_adapter_seq='ATAACAAACGCTAATTGCAA',
                    junk_reads=0, random_reads=0, chimeras=0,
                    glitch_rate=0, glitch_size=0, glitch_skip=0)
        lengths = []
        for _ in range(self.trials):
            fragment, info = \
                badread.simulate.build_fragment(self.lengths, self.ref_seqs, self.rev_comp_ref_seqs,
                                                self.ref_contigs, self.ref_contig_weights,
                                                self.ref_circular, args, start_adapt_rate,
                                                start_adapt_amount, end_adapt_rate,
                                                end_adapt_amount)
            self.assertTrue(fragment.startswith(args.start_adapter_seq))
            self.assertTrue(fragment.endswith(args.end_adapter_seq))
            self.assertLessEqual(len(fragment),
                                 1000 + len(args.start_adapter_seq) + len(args.end_adapter_seq))
            lengths.append(len(fragment))
        most_common_length = max(set(lengths), key=lengths.count)
        self.assertEqual(most_common_length,
                         1000 + len(args.start_adapter_seq) + len(args.end_adapter_seq))

    def test_chimeras(self):
        # High chimera rate.
        start_adapt_rate, start_adapt_amount = 0.0, 0.0
        end_adapt_rate, end_adapt_amount = 0.0, 0.0
        Args = collections.namedtuple('Args', ['start_adapter_seq', 'end_adapter_seq',
                                               'junk_reads', 'random_reads', 'chimeras',
                                               'glitch_rate', 'glitch_size', 'glitch_skip'])
        args = Args(start_adapter_seq='', end_adapter_seq='',
                    junk_reads=0, random_reads=0, chimeras=25,
                    glitch_rate=0, glitch_size=0, glitch_skip=0)
        lengths = []
        for _ in range(self.trials):
            fragment, info = \
                badread.simulate.build_fragment(self.lengths, self.ref_seqs, self.rev_comp_ref_seqs,
                                                self.ref_contigs, self.ref_contig_weights,
                                                self.ref_circular, args, start_adapt_rate,
                                                start_adapt_amount, end_adapt_rate,
                                                end_adapt_amount)
            lengths.append(len(fragment))
        self.assertGreater(max(lengths), 1000)  # Chimeras make for some longer fragments


class TestCircularFragments(unittest.TestCase):

    def setUp(self):
        self.null = open(os.devnull, 'w')
        self.lengths = badread.fragment_lengths.FragmentLengths(1000, 0, self.null)
        self.ref_seqs = {'r': badread.misc.get_random_sequence(10000)}
        self.ref_depths = {'r': 1.0}
        self.ref_circular = {'r': True}
        self.rev_comp_ref_seqs = {name: badread.misc.reverse_complement(seq)
                                  for name, seq in self.ref_seqs.items()}
        self.ref_contigs, self.ref_contig_weights = \
            badread.simulate.get_ref_contig_weights(self.ref_seqs, self.ref_depths)
        self.trials = 100

    def tearDown(self):
        self.null.close()

    def test_plain(self):
        # No adapters, glitches, chimeras.
        start_adapt_rate, start_adapt_amount = 0.0, 0.0
        end_adapt_rate, end_adapt_amount = 0.0, 0.0
        Args = collections.namedtuple('Args', ['start_adapter_seq', 'end_adapter_seq',
                                               'junk_reads', 'random_reads', 'chimeras',
                                               'glitch_rate', 'glitch_size', 'glitch_skip'])
        args = Args(start_adapter_seq='', end_adapter_seq='',
                    junk_reads=0, random_reads=0, chimeras=0,
                    glitch_rate=0, glitch_size=0, glitch_skip=0)
        forward_count, reverse_count = 0, 0
        forward_ref = self.ref_seqs['r'] + self.ref_seqs['r']
        reverse_ref = self.rev_comp_ref_seqs['r'] + self.rev_comp_ref_seqs['r']
        for _ in range(self.trials):
            fragment, info = \
                badread.simulate.build_fragment(self.lengths, self.ref_seqs, self.rev_comp_ref_seqs,
                                                self.ref_contigs, self.ref_contig_weights,
                                                self.ref_circular, args, start_adapt_rate,
                                                start_adapt_amount, end_adapt_rate,
                                                end_adapt_amount)
            if fragment in forward_ref:
                forward_count += 1
            if fragment in reverse_ref:
                reverse_count += 1
            self.assertEqual(len(fragment), 1000)

        # Should have some from each strand.
        self.assertGreater(forward_count, 10)
        self.assertGreater(reverse_count, 10)
        self.assertEqual(forward_count + reverse_count, self.trials)

    def test_full_adapters(self):
        # Start and end adapters reliably put on entirely.
        start_adapt_rate, start_adapt_amount = 1.0, 1.0
        end_adapt_rate, end_adapt_amount = 1.0, 1.0
        Args = collections.namedtuple('Args', ['start_adapter_seq', 'end_adapter_seq',
                                               'junk_reads', 'random_reads', 'chimeras',
                                               'glitch_rate', 'glitch_size', 'glitch_skip'])
        args = Args(start_adapter_seq='TGTATAATACGACGCCGAGC',
                    end_adapter_seq='ATAACAAACGCTAATTGCAA',
                    junk_reads=0, random_reads=0, chimeras=0,
                    glitch_rate=0, glitch_size=0, glitch_skip=0)
        for _ in range(self.trials):
            fragment, info = \
                badread.simulate.build_fragment(self.lengths, self.ref_seqs, self.rev_comp_ref_seqs,
                                                self.ref_contigs, self.ref_contig_weights,
                                                self.ref_circular, args, start_adapt_rate,
                                                start_adapt_amount, end_adapt_rate,
                                                end_adapt_amount)
            self.assertTrue(fragment.startswith(args.start_adapter_seq))
            self.assertTrue(fragment.endswith(args.end_adapter_seq))
            self.assertEqual(len(fragment),
                             1000 + len(args.start_adapter_seq) + len(args.end_adapter_seq))

    def test_chimeras(self):
        # High chimera rate.
        start_adapt_rate, start_adapt_amount = 0.0, 0.0
        end_adapt_rate, end_adapt_amount = 0.0, 0.0
        Args = collections.namedtuple('Args', ['start_adapter_seq', 'end_adapter_seq',
                                               'junk_reads', 'random_reads', 'chimeras',
                                               'glitch_rate', 'glitch_size', 'glitch_skip'])
        args = Args(start_adapter_seq='', end_adapter_seq='',
                    junk_reads=0, random_reads=0, chimeras=25,
                    glitch_rate=0, glitch_size=0, glitch_skip=0)
        lengths = []
        for _ in range(self.trials):
            fragment, info = \
                badread.simulate.build_fragment(self.lengths, self.ref_seqs, self.rev_comp_ref_seqs,
                                                self.ref_contigs, self.ref_contig_weights,
                                                self.ref_circular, args, start_adapt_rate,
                                                start_adapt_amount, end_adapt_rate,
                                                end_adapt_amount)
            lengths.append(len(fragment))
            self.assertTrue(len(fragment) % 1000 == 0)
        self.assertGreater(statistics.mean(lengths), 1000)  # Chimeras make for longer fragments

    def test_enlarging_glitches(self):
        # Glitches with size but no skip make fragments longer.
        start_adapt_rate, start_adapt_amount = 0.0, 0.0
        end_adapt_rate, end_adapt_amount = 0.0, 0.0
        Args = collections.namedtuple('Args', ['start_adapter_seq', 'end_adapter_seq',
                                               'junk_reads', 'random_reads', 'chimeras',
                                               'glitch_rate', 'glitch_size', 'glitch_skip'])
        args = Args(start_adapter_seq='', end_adapter_seq='',
                    junk_reads=0, random_reads=0, chimeras=0,
                    glitch_rate=100, glitch_size=10, glitch_skip=0)
        lengths = []
        for _ in range(self.trials):
            fragment, info = \
                badread.simulate.build_fragment(self.lengths, self.ref_seqs, self.rev_comp_ref_seqs,
                                                self.ref_contigs, self.ref_contig_weights,
                                                self.ref_circular, args, start_adapt_rate,
                                                start_adapt_amount, end_adapt_rate,
                                                end_adapt_amount)
            lengths.append(len(fragment))
        self.assertGreater(statistics.mean(lengths), 1000)

    def test_shrinking_glitches(self):
        # Glitches with skip but no size make fragments shorter.
        start_adapt_rate, start_adapt_amount = 0.0, 0.0
        end_adapt_rate, end_adapt_amount = 0.0, 0.0
        Args = collections.namedtuple('Args', ['start_adapter_seq', 'end_adapter_seq',
                                               'junk_reads', 'random_reads', 'chimeras',
                                               'glitch_rate', 'glitch_size', 'glitch_skip'])
        args = Args(start_adapter_seq='', end_adapter_seq='',
                    junk_reads=0, random_reads=0, chimeras=0,
                    glitch_rate=100, glitch_size=0, glitch_skip=10)
        lengths = []
        for _ in range(self.trials):
            fragment, info = \
                badread.simulate.build_fragment(self.lengths, self.ref_seqs, self.rev_comp_ref_seqs,
                                                self.ref_contigs, self.ref_contig_weights,
                                                self.ref_circular, args, start_adapt_rate,
                                                start_adapt_amount, end_adapt_rate,
                                                end_adapt_amount)
            lengths.append(len(fragment))
        self.assertLess(statistics.mean(lengths), 1000)

    def test_junk_and_random(self):
        start_adapt_rate, start_adapt_amount = 0.0, 0.0
        end_adapt_rate, end_adapt_amount = 0.0, 0.0
        Args = collections.namedtuple('Args', ['start_adapter_seq', 'end_adapter_seq',
                                               'junk_reads', 'random_reads', 'chimeras',
                                               'glitch_rate', 'glitch_size', 'glitch_skip'])
        args = Args(start_adapter_seq='', end_adapter_seq='',
                    junk_reads=25, random_reads=25, chimeras=0,
                    glitch_rate=0, glitch_size=0, glitch_skip=0)
        forward_count, reverse_count = 0, 0
        forward_ref = self.ref_seqs['r'] + self.ref_seqs['r']
        reverse_ref = self.rev_comp_ref_seqs['r'] + self.rev_comp_ref_seqs['r']
        for _ in range(self.trials):
            fragment, info = \
                badread.simulate.build_fragment(self.lengths, self.ref_seqs, self.rev_comp_ref_seqs,
                                                self.ref_contigs, self.ref_contig_weights,
                                                self.ref_circular, args, start_adapt_rate,
                                                start_adapt_amount, end_adapt_rate,
                                                end_adapt_amount)
            if fragment in forward_ref:
                forward_count += 1
            if fragment in reverse_ref:
                reverse_count += 1
            self.assertEqual(len(fragment), 1000)

        # Should have some from each strand.
        self.assertGreater(forward_count, 5)
        self.assertGreater(reverse_count, 5)
        self.assertLess(forward_count + reverse_count, self.trials)


class TestSmallPlasmidBias(unittest.TestCase):

    def setUp(self):
        # Frag lengths of 10000 and reference length 1000.
        self.null = open(os.devnull, 'w')
        self.lengths = badread.fragment_lengths.FragmentLengths(10000, 0, self.null)
        self.ref_seqs = {'r': badread.misc.get_random_sequence(1000)}
        self.ref_depths = {'r': 1.0}
        self.ref_circular = {'r': True}
        self.rev_comp_ref_seqs = {name: badread.misc.reverse_complement(seq)
                                  for name, seq in self.ref_seqs.items()}
        self.ref_contigs, self.ref_contig_weights = \
            badread.simulate.get_ref_contig_weights(self.ref_seqs, self.ref_depths)
        self.trials = 100

    def tearDown(self):
        self.null.close()

    def test_bias_on(self):
        # When bias is on, the program quits with an error, because it can't make any fragments.
        start_adapt_rate, start_adapt_amount = 0.0, 0.0
        end_adapt_rate, end_adapt_amount = 0.0, 0.0
        Args = collections.namedtuple('Args', ['start_adapter_seq', 'end_adapter_seq',
                                               'junk_reads', 'random_reads', 'chimeras',
                                               'glitch_rate', 'glitch_size', 'glitch_skip',
                                               'small_plasmid_bias'])
        args = Args(start_adapter_seq='', end_adapter_seq='',
                    junk_reads=0, random_reads=0, chimeras=0,
                    glitch_rate=0, glitch_size=0, glitch_skip=0,
                    small_plasmid_bias=True)
        with self.assertRaises(SystemExit):
            fragment, info = \
                badread.simulate.build_fragment(self.lengths, self.ref_seqs, self.rev_comp_ref_seqs,
                                                self.ref_contigs, self.ref_contig_weights,
                                                self.ref_circular, args, start_adapt_rate,
                                                start_adapt_amount, end_adapt_rate,
                                                end_adapt_amount)
            self.assertEqual(len(fragment), 1000)


class TestWholeRef(unittest.TestCase):

    def setUp(self):
        # Frag lengths of 10000 and reference length 1000.
        self.null = open(os.devnull, 'w')
        self.lengths = badread.fragment_lengths.FragmentLengths(10000, 0, self.null)
        self.ref_seqs = {'r': badread.misc.get_random_sequence(1000)}
        self.ref_depths = {'r': 1.0}
        self.ref_circular = {'r': False}
        self.rev_comp_ref_seqs = {name: badread.misc.reverse_complement(seq)
                                  for name, seq in self.ref_seqs.items()}
        self.ref_contigs, self.ref_contig_weights = \
            badread.simulate.get_ref_contig_weights(self.ref_seqs, self.ref_depths)
        self.trials = 100

    def tearDown(self):
        self.null.close()

    def test_whole_ref(self):
        # Each fragment should be the whole reference.
        start_adapt_rate, start_adapt_amount = 0.0, 0.0
        end_adapt_rate, end_adapt_amount = 0.0, 0.0
        Args = collections.namedtuple('Args', ['start_adapter_seq', 'end_adapter_seq',
                                               'junk_reads', 'random_reads', 'chimeras',
                                               'glitch_rate', 'glitch_size', 'glitch_skip'])
        args = Args(start_adapter_seq='', end_adapter_seq='',
                    junk_reads=0, random_reads=0, chimeras=0,
                    glitch_rate=0, glitch_size=0, glitch_skip=0)
        forward_count, reverse_count = 0, 0
        for _ in range(self.trials):
            fragment, info = \
                badread.simulate.build_fragment(self.lengths, self.ref_seqs, self.rev_comp_ref_seqs,
                                                self.ref_contigs, self.ref_contig_weights,
                                                self.ref_circular, args, start_adapt_rate,
                                                start_adapt_amount, end_adapt_rate,
                                                end_adapt_amount)
            self.assertEqual(len(fragment), 1000)
            if fragment == self.ref_seqs['r']:
                forward_count += 1
            if fragment == self.rev_comp_ref_seqs['r']:
                reverse_count += 1

        # Should have some from each strand.
        self.assertGreater(forward_count, 10)
        self.assertGreater(reverse_count, 10)
        self.assertEqual(forward_count + reverse_count, self.trials)


class TestRandomJunk(unittest.TestCase):

    def setUp(self):
        self.trials = 20
        self.seq_len = 10000

    def test_random(self):
        # Random sequences compress normally - they get smaller but not too small.
        for _ in range(self.trials):
            random_seq = badread.misc.get_random_sequence(self.seq_len)
            compressed_seq = zlib.compress(random_seq.encode())
            self.assertGreater(len(compressed_seq), len(random_seq) / 10)

    def test_junk(self):
        # Junk sequences compress a lot - they get very small
        for _ in range(self.trials):
            random_seq = badread.simulate.get_junk_fragment(self.seq_len)
            compressed_seq = zlib.compress(random_seq.encode())
            self.assertLess(len(compressed_seq), len(random_seq) / 10)
