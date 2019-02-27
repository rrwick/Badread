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

import os
import unittest

import badread.alignment


class TestAlignments(unittest.TestCase):

    def setUp(self):
        self.ref = os.path.join(os.path.dirname(__file__), 'test_alignment_ref.fasta')
        self.reads = os.path.join(os.path.dirname(__file__), 'test_alignment_reads.fasta')
        self.paf = os.path.join(os.path.dirname(__file__), 'test_alignment.paf')
        self.alignments = []
        with open(self.paf) as paf:
            for line in paf:
                self.alignments.append(badread.alignment.Alignment(line))
        self.null = open(os.devnull, 'w')

    def tearDown(self):
        self.null.close()

    def test_read_names(self):
        self.assertEqual(self.alignments[0].read_name, 'read_1')
        self.assertEqual(self.alignments[1].read_name, 'read_2')

    def test_repr(self):
        self.assertEqual(str(self.alignments[0]), 'read_1:0-1438(+),ref:754-2191(99.861%)')
        self.assertEqual(str(self.alignments[1]), 'read_2:0-1128(-),ref:3357-4486(99.911%)')

    def test_identity(self):
        self.assertGreater(self.alignments[0].percent_identity, 99)
        self.assertGreater(self.alignments[1].percent_identity, 99)

    def test_bad_paf_1(self):
        with self.assertRaises(SystemExit):
            badread.alignment.Alignment('this_is_not_a_paf_line')

    def test_bad_paf_2(self):
        # Missing alignment score
        with self.assertRaises(SystemExit):
            badread.alignment.Alignment('read_1\t1438\t0\t1438\t+\tref\t10000\t754\t2191\t1436\t'
                                        '1438\t60\tNM:i:2\tms:i:2862\tnn:i:0\ttp:A:P\t'
                                        'cm:i:254\ts1:i:1410\ts2:i:0\tdv:f:0.0016\tcg:Z:798M1I639M')

    def test_bad_paf_3(self):
        # Missing CIGAR
        with self.assertRaises(SystemExit):
            badread.alignment.Alignment('read_1\t1438\t0\t1438\t+\tref\t10000\t754\t2191\t1436\t'
                                        '1438\t60\tNM:i:2\tms:i:2862\tAS:i:2862\tnn:i:0\ttp:A:P\t'
                                        'cm:i:254\ts1:i:1410\ts2:i:0\tdv:f:0.0016')


class TestLoadAlignments(unittest.TestCase):

    def setUp(self):
        self.paf = os.path.join(os.path.dirname(__file__), 'test_alignment.paf')
        self.null = open(os.devnull, 'w')

    def tearDown(self):
        self.null.close()

    def test_load(self):
        alignments = badread.alignment.load_alignments(self.paf, output=self.null, dot_interval=1)
        self.assertEqual(len(alignments), 2)

    def test_load_with_max(self):
        alignments = badread.alignment.load_alignments(self.paf, output=self.null, max_alignments=1)
        self.assertEqual(len(alignments), 1)


class TestAlignSequences(unittest.TestCase):

    def test_align_sequences_1(self):
        unaligned_read = 'ACTACGCACTACG'
        unaligned_qual = '1253176253763'
        unaligned_ref = 'ACTACGACTACG'
        alignment = badread.alignment.Alignment('read_1\t13\t0\t13\t+\tref\t12\t0\t12\t12\t13\t255'
                                                '\tAS:i:2862\tcg:Z:6M1I6M')
        aligned_read_seq, aligned_read_qual, aligned_ref_seq, errors_per_read_pos = \
            badread.alignment.align_sequences(unaligned_read, unaligned_qual, unaligned_ref,
                                              alignment)
        self.assertEqual(aligned_read_seq, 'ACTACGCACTACG')
        self.assertEqual(aligned_read_qual, '1253176253763')
        self.assertEqual(aligned_ref_seq, 'ACTACG-ACTACG')
        self.assertEqual(errors_per_read_pos, [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0])

    def test_align_sequences_2(self):
        unaligned_read = 'ACTACGACAACG'
        unaligned_qual = '125317253763'
        unaligned_ref = 'ACTACGCACTACG'
        alignment = badread.alignment.Alignment('read_1\t12\t0\t12\t+\tref\t13\t0\t13\t11\t13\t255'
                                                '\tAS:i:2862\tcg:Z:6M1D6M')
        aligned_read_seq, aligned_read_qual, aligned_ref_seq, errors_per_read_pos = \
            badread.alignment.align_sequences(unaligned_read, unaligned_qual, unaligned_ref,
                                              alignment, gap_char=' ')
        self.assertEqual(aligned_read_seq, 'ACTACG ACAACG')
        self.assertEqual(aligned_read_qual, '125317 253763')
        self.assertEqual(aligned_ref_seq, 'ACTACGCACTACG')
        self.assertEqual(errors_per_read_pos, [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0])
