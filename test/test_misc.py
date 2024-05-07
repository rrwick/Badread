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

import gzip
import os
import unittest

import badread.misc


class TestCompressionType(unittest.TestCase):

    def test_no_compression(self):
        filename = os.path.join(os.path.dirname(__file__), 'test_ref_1.fasta')
        comp_type = badread.misc.get_compression_type(filename)
        self.assertEqual(comp_type, 'plain')

    def test_gzip(self):
        filename = os.path.join(os.path.dirname(__file__), 'test_ref_1.fasta.gz')
        comp_type = badread.misc.get_compression_type(filename)
        self.assertEqual(comp_type, 'gz')

    def test_bzip2(self):
        filename = os.path.join(os.path.dirname(__file__), 'test_ref_1.fasta.bz2')
        with self.assertRaises(SystemExit) as cm:
            badread.misc.get_compression_type(filename)
        self.assertTrue('cannot use bzip2 format' in str(cm.exception))

    def test_zip(self):
        filename = os.path.join(os.path.dirname(__file__), 'test_ref_1.fasta.zip')
        with self.assertRaises(SystemExit) as cm:
            badread.misc.get_compression_type(filename)
        self.assertTrue('cannot use zip format' in str(cm.exception))

    def test_regular_open(self):
        filename = os.path.join(os.path.dirname(__file__), 'test_ref_1.fasta')
        open_func = badread.misc.get_open_func(filename)
        self.assertEqual(open_func, open)

    def test_gzip_open(self):
        filename = os.path.join(os.path.dirname(__file__), 'test_ref_1.fasta.gz')
        open_func = badread.misc.get_open_func(filename)
        self.assertEqual(open_func, gzip.open)


class TestSequenceFileType(unittest.TestCase):

    def test_fastq(self):
        filename = os.path.join(os.path.dirname(__file__), 'test_reads_1.fastq')
        self.assertEqual(badread.misc.get_sequence_file_type(filename), 'FASTQ')

    def test_fasta(self):
        filename = os.path.join(os.path.dirname(__file__), 'test_ref_1.fasta')
        self.assertEqual(badread.misc.get_sequence_file_type(filename), 'FASTA')

    def test_fasta_gz(self):
        filename = os.path.join(os.path.dirname(__file__), 'test_ref_1.fasta.gz')
        self.assertEqual(badread.misc.get_sequence_file_type(filename), 'FASTA')

    def test_neither(self):
        filename = os.path.join(os.path.dirname(__file__), 'test_alignment.paf')
        with self.assertRaises(ValueError):
            badread.misc.get_sequence_file_type(filename)

    def test_not_a_file(self):
        filename = os.path.join(os.path.dirname(__file__), 'sdjkhfksdjhfksjd')
        with self.assertRaises(SystemExit):
            badread.misc.get_sequence_file_type(filename)

    def test_empty_file(self):
        filename = os.path.join(os.path.dirname(__file__), 'test_empty_file')
        with self.assertRaises(ValueError):
            badread.misc.get_sequence_file_type(filename)

    def test_bad_file(self):
        filename = os.path.join(os.path.dirname(__file__), 'test_bad_file')
        with self.assertRaises(ValueError):
            badread.misc.get_sequence_file_type(filename)


class TestReverseComplement(unittest.TestCase):

    def test_comp_base_1(self):
        self.assertEqual(badread.misc.complement_base('A'), 'T')

    def test_comp_base_2(self):
        self.assertEqual(badread.misc.complement_base('C'), 'G')

    def test_comp_base_3(self):
        self.assertEqual(badread.misc.complement_base('G'), 'C')

    def test_comp_base_4(self):
        self.assertEqual(badread.misc.complement_base('T'), 'A')

    def test_comp_base_5(self):
        self.assertEqual(badread.misc.complement_base('N'), 'N')

    def test_comp_base_6(self):
        self.assertEqual(badread.misc.complement_base('%'), 'N')

    def test_rev_comp_1(self):
        self.assertEqual(badread.misc.reverse_complement('AATTCC'), 'GGAATT')

    def test_rev_comp_2(self):
        self.assertEqual(badread.misc.reverse_complement('gAgA'), 'TcTc')


class TestLoadSequences(unittest.TestCase):

    def setUp(self):
        self.null = open(os.devnull, 'w')

    def tearDown(self):
        self.null.close()

    def check_fasta(self, seqs, depths, circular):
        self.assertEqual(len(seqs), 11)
        self.assertEqual(seqs['A'], 'GTCCGCTCAACCTCCTGTTTATCTTAAACC')
        self.assertEqual(depths['A'], 1.0)
        self.assertEqual(circular['A'], True)
        self.assertEqual(seqs['E'], 'ATAAGTCCCGCCCTTCGCCCGTTTTCAGGG')
        self.assertEqual(depths['E'], 3.0)
        self.assertEqual(circular['E'], False)

    def check_fastq(self, reads):
        self.assertEqual(len(reads), 3)
        self.assertTrue(reads['read_1'][0].startswith('GTTACTTCGATATC'))
        self.assertTrue(reads['read_1'][1].startswith('BCCCBGGGGGGCGG'))
        self.assertTrue(reads['read_2'][0].startswith('CACATACAGGCAGA'))
        self.assertTrue(reads['read_2'][1].startswith('CCCCCFG>GG1GGG'))
        self.assertTrue(reads['read_3'][0].startswith('ACAAATATCAGGAT'))
        self.assertTrue(reads['read_3'][1].startswith('BABCBGGGGGGGGG'))

    def test_load_non_gzipped_fasta(self):
        filename = os.path.join(os.path.dirname(__file__), 'test_ref_1.fasta')
        seqs, depths, circular = badread.misc.load_fasta(filename)
        self.check_fasta(seqs, depths, circular)

    def test_load_gzipped_fasta(self):
        filename = os.path.join(os.path.dirname(__file__), 'test_ref_1.fasta.gz')
        seqs, depths, circular = badread.misc.load_fasta(filename)
        self.check_fasta(seqs, depths, circular)

    def test_load_bad_format_fasta(self):
        filename = os.path.join(os.path.dirname(__file__), 'test_ref_1_bad.fasta')
        seqs, depths, circular = badread.misc.load_fasta(filename)
        self.check_fasta(seqs, depths, circular)

    def test_load_fastq(self):
        filename = os.path.join(os.path.dirname(__file__), 'test_reads_1.fastq')
        reads = badread.misc.load_fastq(filename, output=self.null)
        self.check_fastq(reads)

    def test_load_bad_format_fastq(self):
        filename = os.path.join(os.path.dirname(__file__), 'test_reads_2.fastq')
        reads = badread.misc.load_fastq(filename, output=self.null, dot_interval=1)
        self.check_fastq(reads)

    def test_load_fastq_wrong_type(self):
        filename = os.path.join(os.path.dirname(__file__), 'test_ref_1.fasta')
        with self.assertRaises(SystemExit) as cm:
            badread.misc.load_fastq(filename, output=self.null)
        self.assertTrue('is not FASTQ format' in str(cm.exception))


class TestRandomSeqs(unittest.TestCase):

    def test_random_base(self):
        for _ in range(10):
            b = badread.misc.get_random_base()
            self.assertTrue(b == 'A' or b == 'C' or b == 'G' or b == 'T')

    def test_random_different_base(self):
        for b in ['A', 'C', 'G', 'T']:
            for _ in range(10):
                diff_b = badread.misc.get_random_different_base(b)
                self.assertTrue(b == 'A' or b == 'C' or b == 'G' or b == 'T')
                self.assertNotEqual(b, diff_b)

    def test_random_seq(self):
        for seq_len in range(100):
            random_seq = badread.misc.get_random_sequence(seq_len)
            self.assertEqual(len(random_seq), seq_len)


class TestNumFormatting(unittest.TestCase):

    def test_float_to_str_1(self):
        self.assertEqual(badread.misc.float_to_str(1.00), '1')

    def test_float_to_str_2(self):
        self.assertEqual(badread.misc.float_to_str(1.01, decimals=2), '1.01')

    def test_float_to_str_3(self):
        self.assertEqual(badread.misc.float_to_str(1.01, decimals=1), '1.0')

    def test_float_to_str_4(self):
        self.assertEqual(badread.misc.float_to_str(1.06, decimals=1), '1.1')

    def test_float_to_str_5(self):
        self.assertEqual(badread.misc.float_to_str(1.01, decimals=4), '1.0100')

    def test_float_to_str_6(self):
        self.assertEqual(badread.misc.float_to_str(1.01, decimals=4, trim_zeros=True), '1.01')


class TestOther(unittest.TestCase):

    def test_str_is_int_1(self):
        self.assertTrue(badread.misc.str_is_int('12'))

    def test_str_is_int_2(self):
        self.assertFalse(badread.misc.str_is_int('abc'))

    def test_str_is_dna_sequence_1(self):
        self.assertTrue(badread.misc.str_is_dna_sequence('ACGACTCAGACT'))

    def test_str_is_dna_sequence_2(self):
        self.assertFalse(badread.misc.str_is_dna_sequence('abcdefg'))

    def test_identity_from_edlib_cigar_1(self):
        self.assertEqual(badread.misc.identity_from_edlib_cigar('10='), 1.0)

    def test_identity_from_edlib_cigar_2(self):
        self.assertEqual(badread.misc.identity_from_edlib_cigar('10D'), 0.0)

    def test_identity_from_edlib_cigar_3(self):
        self.assertEqual(badread.misc.identity_from_edlib_cigar('5=5X'), 0.5)

    def test_identity_from_edlib_cigar_4(self):
        self.assertEqual(badread.misc.identity_from_edlib_cigar(''), 0.0)

    def test_random_chance_0(self):
        successes = sum(1 if badread.misc.random_chance(0.0) else 0 for _ in range(1000))
        self.assertEqual(successes, 0)

    def test_random_chance_100(self):
        successes = sum(1 if badread.misc.random_chance(1.0) else 0 for _ in range(1000))
        self.assertEqual(successes, 1000)

    def test_random_chance_50(self):
        successes = sum(1 if badread.misc.random_chance(0.5) else 0 for _ in range(1000))
        self.assertTrue(200 < successes < 800)

    def test_only_acgt(self):
        self.assertTrue(badread.misc.only_acgt("AACGATCAGCACTG"))
        self.assertTrue(badread.misc.only_acgt("CGCGCGCGCGCG"))
        self.assertTrue(badread.misc.only_acgt("TTTTTTTT"))
        self.assertTrue(badread.misc.only_acgt("A"))
        self.assertTrue(badread.misc.only_acgt(""))
        self.assertFalse(badread.misc.only_acgt("ACGANCTCG"))
        self.assertFalse(badread.misc.only_acgt("CGCTXACGACT"))
        self.assertFalse(badread.misc.only_acgt("12345"))
        self.assertFalse(badread.misc.only_acgt("acgactacgac"))
        self.assertFalse(badread.misc.only_acgt("ACGAtCGACG"))
