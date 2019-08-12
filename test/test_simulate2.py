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

from io import StringIO
import collections
import os
import unittest

import badread.simulate
import badread.identities
import badread.error_model
import badread.qscore_model
import badread.misc


def sequence(reference_filename, read_count=5000, mean_frag_length=100, small_plasmid_bias=False,
             seed=None, mean_identity=85):
    quantity = mean_frag_length * read_count
    Args = collections.namedtuple('Args', ['reference', 'quantity',
                                           'mean_frag_length', 'frag_length_stdev',
                                           'mean_identity', 'max_identity', 'identity_stdev',
                                           'error_model', 'qscore_model', 'seed',
                                           'start_adapter', 'end_adapter',
                                           'start_adapter_seq', 'end_adapter_seq',
                                           'junk_reads', 'random_reads', 'chimeras',
                                           'glitch_rate', 'glitch_size', 'glitch_skip',
                                           'small_plasmid_bias'])
    args = Args(reference=reference_filename, quantity=quantity,
                mean_frag_length=mean_frag_length, frag_length_stdev=10,
                mean_identity=mean_identity, max_identity=95, identity_stdev=5,
                error_model='random', qscore_model='ideal', seed=seed,
                start_adapter='0,0', end_adapter='0,0',
                start_adapter_seq='', end_adapter_seq='',
                junk_reads=0, random_reads=0, chimeras=0,
                glitch_rate=0, glitch_size=0, glitch_skip=0,
                small_plasmid_bias=small_plasmid_bias)

    with open(os.devnull, 'w') as null:
        badread.simulate.simulate(args, output=null)


def count_strands(reads):
    forward_count, reverse_count = 0, 0
    s = StringIO(reads)
    for line in s:
        header = line.strip()
        if ',+strand,' in header:
            forward_count += 1
        elif ',-strand,' in header:
            reverse_count += 1
        else:
            assert False
        next(s)
        next(s)
        next(s)
    return forward_count, reverse_count


def count_ref_contigs(reads):
    counts = collections.defaultdict(int)
    s = StringIO(reads)
    for line in s:
        header = line.strip()
        name = header.split(' ')[1].split(',')[0]
        counts[name] += 1
        next(s)
        next(s)
        next(s)
    return counts


class TestSimulate(unittest.TestCase):

    def test_strand_count(self):
        # Should get reads from both strands.
        ref_filename = os.path.join(os.path.dirname(__file__), 'test_ref_2.fasta')
        with badread.misc.captured_output() as (out, err):
            sequence(ref_filename)
        out, err = out.getvalue().strip(), err.getvalue().strip()
        forward_count, reverse_count = count_strands(out)
        total = forward_count + reverse_count
        self.assertGreater(forward_count, total / 4)
        self.assertGreater(reverse_count, total / 4)

    def test_counts_1(self):
        # Contigs A and B have the same length, but B is double the depth so should get 2/3 of
        # the reads.
        ref_filename = os.path.join(os.path.dirname(__file__), 'test_ref_2.fasta')
        with badread.misc.captured_output() as (out, err):
            sequence(ref_filename)
        out, err = out.getvalue().strip(), err.getvalue().strip()
        counts = count_ref_contigs(out)
        a_count, b_count = counts['A'], counts['B']
        ratio = b_count / a_count
        self.assertGreater(ratio, 1.6)
        self.assertLess(ratio, 2.4)

    def test_counts_2(self):
        # Contigs C and D have the same depth, but C is double the length so should get 2/3 of
        # the reads.
        ref_filename = os.path.join(os.path.dirname(__file__), 'test_ref_3.fasta')
        with badread.misc.captured_output() as (out, err):
            sequence(ref_filename)
        out, err = out.getvalue().strip(), err.getvalue().strip()
        counts = count_ref_contigs(out)
        c_count, d_count = counts['C'], counts['D']
        ratio = c_count / d_count
        self.assertGreater(ratio, 1.6)
        self.assertLess(ratio, 2.4)

    def test_counts_3(self):
        # Contig E is half the length but 4x the depth of contig F, so it should get 2/3 of the
        # reads.
        ref_filename = os.path.join(os.path.dirname(__file__), 'test_ref_4.fasta')
        with badread.misc.captured_output() as (out, err):
            sequence(ref_filename)
        out, err = out.getvalue().strip(), err.getvalue().strip()
        counts = count_ref_contigs(out)
        e_count, f_count = counts['E'], counts['F']
        ratio = e_count / f_count
        self.assertGreater(ratio, 1.6)
        self.assertLess(ratio, 2.4)

    def test_plasmid_bias_on(self):
        # When --small_plasmid_bias is on, contig H gets few to no reads, because its length
        # (50 bp) is much less than the mean fragment length (100 bp).
        ref_filename = os.path.join(os.path.dirname(__file__), 'test_ref_5.fasta')
        with badread.misc.captured_output() as (out, err):
            sequence(ref_filename, small_plasmid_bias=True)
        out, err = out.getvalue().strip(), err.getvalue().strip()
        counts = count_ref_contigs(out)
        g_count, h_count = counts['G'], counts['H']
        ratio = h_count / g_count
        self.assertLess(ratio, 0.2)

    def test_no_seed(self):
        # Without a seed, repeated runs should give different results.
        ref_filename = os.path.join(os.path.dirname(__file__), 'test_ref_2.fasta')
        with badread.misc.captured_output() as (out1, err1):
            sequence(ref_filename, read_count=10, seed=None)
        out1, err1 = out1.getvalue().strip(), err1.getvalue().strip()
        with badread.misc.captured_output() as (out2, err2):
            sequence(ref_filename, read_count=10, seed=None)
        out2, err2 = out2.getvalue().strip(), err2.getvalue().strip()
        self.assertNotEqual(out1, out2)

    def test_with_seed(self):
        # With a seed, repeated runs should give the same results.
        ref_filename = os.path.join(os.path.dirname(__file__), 'test_ref_2.fasta')
        with badread.misc.captured_output() as (out1, err1):
            sequence(ref_filename, read_count=10, seed=1)
        out1, err1 = out1.getvalue().strip(), err1.getvalue().strip()
        with badread.misc.captured_output() as (out2, err2):
            sequence(ref_filename, read_count=10, seed=1)
        out2, err2 = out2.getvalue().strip(), err2.getvalue().strip()
        with badread.misc.captured_output() as (out3, err3):
            sequence(ref_filename, read_count=10, seed=2)
        out3, err3 = out3.getvalue().strip(), err3.getvalue().strip()
        with badread.misc.captured_output() as (out4, err4):
            sequence(ref_filename, read_count=10, seed=2)
        out4, err4 = out4.getvalue().strip(), err4.getvalue().strip()
        self.assertEqual(out1, out2)
        self.assertEqual(out3, out4)
        self.assertNotEqual(out1, out3)

    def test_very_low_id(self):
        # If we ask for reads with extremely low id, Badread should do the best it can (not get
        # caught in an infinite loop).
        ref_filename = os.path.join(os.path.dirname(__file__), 'test_ref_2.fasta')
        with badread.misc.captured_output() as (out, err):
            sequence(ref_filename, mean_identity=30, read_count=100)
        out, err = out.getvalue().strip(), err.getvalue().strip()
        line_count = len(out.splitlines())
        self.assertEqual(line_count % 4, 0)
        self.assertGreater(line_count, 20)
