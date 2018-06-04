"""
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Badread

This module contains some tests for Badread. To run them, execute `python3 -m unittest` from the
root Badread directory.

This file is part of Badread. Badread is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Badread is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Badread.
If not, see <http://www.gnu.org/licenses/>.
"""

import itertools
import os
import unittest
import badread.error_model


class TestKmerAlignment(unittest.TestCase):
    """
    Tests the align_kmers function (k-mer vs alternative alignment).
    """
    def test_equal_1(self):
        self.assertEqual(badread.error_model.align_kmers('ACGT', 'ACGT'),
                         ['A', 'C', 'G', 'T'])

    def test_equal_2(self):
        self.assertEqual(badread.error_model.align_kmers('ACGACTAGCTACG', 'ACGACTAGCTACG'),
                         ['A', 'C', 'G', 'A', 'C', 'T', 'A', 'G', 'C', 'T', 'A', 'C', 'G'])

    def test_substitutions_1(self):
        self.assertEqual(badread.error_model.align_kmers('ACGT', 'ACTT'),
                         ['A', 'C', 'T', 'T'])

    def test_substitutions_2(self):
        self.assertEqual(badread.error_model.align_kmers('CGTC', 'CCTC'),
                         ['C', 'C', 'T', 'C'])

    def test_substitutions_3(self):
        self.assertEqual(badread.error_model.align_kmers('ACGACTAGCTACG', 'ACGACAAGCTACG'),
                         ['A', 'C', 'G', 'A', 'C', 'A', 'A', 'G', 'C', 'T', 'A', 'C', 'G'])

    def test_deletion_1(self):
        self.assertEqual(badread.error_model.align_kmers('ACGT', 'AGT'),
                         ['A', '', 'G', 'T'])

    def test_deletion_2(self):
        self.assertEqual(badread.error_model.align_kmers('ACGT', 'ACT'),
                         ['A', 'C', '', 'T'])

    def test_deletion_3(self):
        self.assertEqual(badread.error_model.align_kmers('ACGACTAGCTACG', 'ACGACTGCTACG'),
                         ['A', 'C', 'G', 'A', 'C', 'T', '', 'G', 'C', 'T', 'A', 'C', 'G'])

    def test_deletion_4(self):
        # When there are multiple possible alignments, all are acceptable.
        result = badread.error_model.align_kmers('AAAA', 'AAA')
        self.assertTrue(result == ['A', 'A', '', 'A'] or
                        result == ['A', '', 'A', 'A'])

    def test_deletion_5(self):
        # When there are multiple possible alignments, all are acceptable.
        result = badread.error_model.align_kmers('CCCCC', 'CCCC')
        self.assertTrue(result == ['C', '', 'C', 'C', 'C'] or
                        result == ['C', 'C', '', 'C', 'C'] or
                        result == ['C', 'C', 'C', '', 'C'])

    def test_deletion_6(self):
        # When there are multiple possible alignments, all are acceptable.
        result = badread.error_model.align_kmers('CCCCC', 'CCC')
        self.assertTrue(result == ['C', '', '', 'C', 'C'] or
                        result == ['C', '', 'C', '', 'C'] or
                        result == ['C', 'C', '', '', 'C'])

    def test_deletion_7(self):
        # The first and last base should always be the same.
        self.assertEqual(badread.error_model.align_kmers('ACGACTAGCTACG', 'AG'),
                         ['A', '', '', '', '', '', '', '', '', '', '', '', 'G'])

    def test_deletion_8(self):
        # The first and last base should always match.
        self.assertEqual(badread.error_model.align_kmers('AAAAAAAAAAAAA', 'AA'),
                         ['A', '', '', '', '', '', '', '', '', '', '', '', 'A'])

    def test_insertion_1(self):
        # Insertions can go before or after a base.
        result = badread.error_model.align_kmers('ACGT', 'ACAGT')
        self.assertTrue(result == ['A', 'CA', 'G', 'T'] or
                        result == ['A', 'C', 'AG', 'T'])

    def test_insertion_2(self):
        # The first and last base should always be unaltered, so if the insertions is at the
        # beginning, there's only one place for it.
        self.assertEqual(badread.error_model.align_kmers('ACGT', 'ATCGT'),
                         ['A', 'TC', 'G', 'T'])

    def test_insertion_3(self):
        # The first and last base should always be unaltered, so if the insertions is at the
        # end, there's only one place for it.
        self.assertEqual(badread.error_model.align_kmers('ACGT', 'ACGCT'),
                         ['A', 'C', 'GC', 'T'])

    def test_insertion_4(self):
        # When there are multiple possible alignments, all are acceptable.
        result = badread.error_model.align_kmers('AAAA', 'AAAAA')
        self.assertTrue(result == ['A', 'A', 'AA', 'A'] or
                        result == ['A', 'AA', 'A', 'A'])

    def test_insertion_5(self):
        # When there are multiple possible alignments, all are acceptable.
        result = badread.error_model.align_kmers('AAAA', 'AAAAAA')
        self.assertTrue(result == ['A', 'A', 'AAA', 'A'] or
                        result == ['A', 'AAA', 'A', 'A'] or
                        result == ['A', 'AA', 'AA', 'A'])

    def test_insertion_6(self):
        # When there are multiple possible alignments, all are acceptable.
        result = badread.error_model.align_kmers('ACGGAT', 'ACGGGAT')
        self.assertTrue(result == ['A', 'CG', 'G', 'G', 'A', 'T'] or
                        result == ['A', 'C', 'GG', 'G', 'A', 'T'] or
                        result == ['A', 'C', 'G', 'GG', 'A', 'T'] or
                        result == ['A', 'C', 'G', 'G', 'GA', 'T'])

    def test_insertion_7(self):
        # The first and last base should always match.
        self.assertEqual(badread.error_model.align_kmers('AAA', 'AAAAAAAAAAAAAA'),
                         ['A', 'AAAAAAAAAAAA', 'A'])

    def test_assertion_fail_1(self):
        # Reference k-mer must be at least a 3-mer.
        with self.assertRaises(AssertionError):
            badread.error_model.align_kmers('AC', 'AC')

    def test_assertion_fail_2(self):
        # Alt k-mer must be at least a 2-mer.
        with self.assertRaises(AssertionError):
            badread.error_model.align_kmers('AAA', 'A')

    def test_assertion_fail_3(self):
        # The first and last base should always match.
        with self.assertRaises(AssertionError):
            badread.error_model.align_kmers('ACGTA', 'ACGTC')


class TestLoadErrorModel(unittest.TestCase):
    """
    Loads a simple 4-mer error model from a file and makes sure it looks okay.
    """
    def setUp(self):
        model_filename = os.path.join(os.path.dirname(__file__), '4-mer_model')
        self.model = badread.error_model.ErrorModel(model_filename)

    def test_k_size(self):
        self.assertEqual(self.model.kmer_size, 4)

    def test_kmers(self):
        self.assertEqual(len(self.model.alternatives), 256)
        self.assertEqual(len(self.model.probabilities), 256)
        all_4_mers = [''.join(x) for x in itertools.product('ACGT', repeat=4)]
        for kmer in all_4_mers:
            self.assertTrue(kmer in self.model.alternatives)
            self.assertTrue(kmer in self.model.probabilities)

    def test_specifics_1(self):
        alts = self.model.alternatives['AAAA']
        probs = self.model.probabilities['AAAA']
        self.assertEqual(len(alts), 6)
        self.assertEqual(len(probs), 6)
        self.assertEqual(alts[5], ['A', 'G', 'A', 'A'])
        self.assertAlmostEqual(probs[5], 0.003496)

    def test_specifics_2(self):
        alts = self.model.alternatives['GCCA']
        probs = self.model.probabilities['GCCA']
        self.assertEqual(len(alts), 6)
        self.assertEqual(len(probs), 6)
        self.assertEqual(alts[5], ['G', 'C', 'CA', 'A'])
        self.assertAlmostEqual(probs[5], 0.005422)


class Test4MerErrorModel(unittest.TestCase):
    """
    Uses a simple 4-mer error model to make errors.
    """
    def setUp(self):
        model_filename = os.path.join(os.path.dirname(__file__), '4-mer_model')
        self.model = badread.error_model.ErrorModel(model_filename)

    def test_ACAC(self):
        # The model never gets this k-mer wrong.
        new_kmer = self.model.add_errors_to_kmer('ACAC')
        for i in range(100):
            self.assertEqual(new_kmer, ['A', 'C', 'A', 'C'])

    def test_ACAG(self):
        # The gets this k-mer wrong half of the time, always to ACGG.
        correct_count, alt_count = 0, 0
        for i in range(10000):
            new_kmer = self.model.add_errors_to_kmer('ACAG')
            self.assertTrue(new_kmer == ['A', 'C', 'A', 'G'] or
                            new_kmer == ['A', 'C', 'G', 'G'])
            if new_kmer == ['A', 'C', 'A', 'G']:
                correct_count += 1
            else:
                alt_count += 1
        self.assertTrue(4000 < correct_count < 6000)
        self.assertTrue(4000 < alt_count < 6000)

    def test_ACAT(self):
        # The gets this k-mer wrong half of the time. Half of those errors are to ACGT and the
        # other half are to ATAT.
        correct_count, alt_1_count, alt_2_count = 0, 0, 0
        for i in range(10000):
            new_kmer = self.model.add_errors_to_kmer('ACAT')
            self.assertTrue(new_kmer == ['A', 'C', 'A', 'T'] or
                            new_kmer == ['A', 'C', 'G', 'T'] or
                            new_kmer == ['A', 'T', 'A', 'T'])
            if new_kmer == ['A', 'C', 'A', 'T']:
                correct_count += 1
            elif new_kmer == ['A', 'C', 'G', 'T']:
                alt_1_count += 1
            elif new_kmer == ['A', 'T', 'A', 'T']:
                alt_2_count += 1
        self.assertTrue(4000 < correct_count < 6000)
        self.assertTrue(2000 < alt_1_count < 3000)
        self.assertTrue(2000 < alt_2_count < 3000)

    def test_ACCA(self):
        # The gets this k-mer wrong half of the time. Half of those errors are to AGGA and the
        # other half are random.
        correct_count, alt_count = 0, 0
        new_kmers = set()
        for i in range(10000):
            new_kmer = self.model.add_errors_to_kmer('ACCA')
            if new_kmer == ['A', 'C', 'C', 'A']:
                correct_count += 1
            elif new_kmer == ['A', 'G', 'G', 'A']:
                alt_count += 1
            new_kmers.add('_'.join(new_kmer))
        self.assertTrue(4000 < correct_count < 6000)
        self.assertTrue(2000 < alt_count < 3000)

        # There are 46 new k-mers in total: 44 from the random one-change k-mers, one from AGGA
        # (which has two changes) and one from the unmutated original.
        self.assertEqual(len(new_kmers), 46)

    def test_ACCC(self):
        # The model always gets this k-mer wrong (as ACC).
        new_kmer = self.model.add_errors_to_kmer('ACCC')
        for i in range(100):
            self.assertEqual(''.join(new_kmer), 'ACC')

    def test_ACCG(self):
        # The model always gets this k-mer wrong: half the time to ACG and half to a random error.
        alt_count = 0
        new_kmers = set()
        for i in range(10000):
            new_kmer = self.model.add_errors_to_kmer('ACCG')
            self.assertNotEqual(new_kmer, ['A', 'C', 'C', 'G'])
            if ''.join(new_kmer) == 'ACG':
                alt_count += 1
            new_kmers.add('_'.join(new_kmer))
        self.assertTrue(4000 < alt_count < 6000)

        # There are 44 new k-mers in total: the random one-change k-mers (ACG is a one-change k-mer
        # and is included, and the correct k-mer never occurs).
        self.assertEqual(len(new_kmers), 44)

    def test_ACCT(self):
        # The model always gets this k-mer correct half the time and a random error half the time.
        correct_count = 0
        new_kmers = set()
        for i in range(10000):
            new_kmer = self.model.add_errors_to_kmer('ACCT')
            if new_kmer == ['A', 'C', 'C', 'T']:
                correct_count += 1
            new_kmers.add('_'.join(new_kmer))
        self.assertTrue(4000 < correct_count < 6000)

        # There are 45 new k-mers in total: 44 random one-change k-mers and the correct k-mer.
        self.assertEqual(len(new_kmers), 45)


class TestRandomErrorModel(unittest.TestCase):
    """
    Tests a random error model (i.e. an error model not based on k-mers and loaded from a file).
    """
    def setUp(self):
        self.model = badread.error_model.ErrorModel('random')

    def test_k_size(self):
        self.assertEqual(self.model.kmer_size, 1)

    def test_1_mer(self):
        """
        There are 11 possible mutated versions of a 1-mer: 3 subs, 1 deletion and 7 insertions
        (actually 8 insertions, but two are the same because when inserting the same base it
        doesn't matter if it goes on the front or the back).
        """
        new_kmers = set()
        for i in range(10000):  # try many times to make sure we get them all
            new_kmer = self.model.add_errors_to_kmer('A')
            new_kmers.add('_'.join(new_kmer))
        self.assertEqual(len(new_kmers), 11)

    def test_2_mer(self):
        """
        There are 16 possible mutated versions of a 2-mer: 6 subs, 2 deletions and 14 insertions.
        """
        new_kmers = set()
        for i in range(10000):  # try many times to make sure we get them all
            new_kmer = self.model.add_errors_to_kmer('AC')
            new_kmers.add('_'.join(new_kmer))
        self.assertEqual(len(new_kmers), 22)

    def test_3_mer(self):
        """
        There are 24 possible mutated versions of a 3-mer: 9 subs, 3 deletions and 21 insertions.
        """
        new_kmers = set()
        for i in range(10000):  # try many times to make sure we get them all
            new_kmer = self.model.add_errors_to_kmer('ACG')
            new_kmers.add('_'.join(new_kmer))
        self.assertEqual(len(new_kmers), 33)

    def test_4_mer(self):
        """
        There are 32 possible mutated versions of a 4-mer: 12 subs, 4 deletions and 28 insertions.
        """
        new_kmers = set()
        for i in range(10000):  # try many times to make sure we get them all
            new_kmer = self.model.add_errors_to_kmer('ACGT')
            new_kmers.add('_'.join(new_kmer))
        self.assertEqual(len(new_kmers), 44)


class TestPerfectErrorModel(unittest.TestCase):
    """
    Tests a perfect error model (no errors ever added to k-mers).
    """
    def setUp(self):
        self.model = badread.error_model.ErrorModel('perfect')

    def test_k_size(self):
        self.assertEqual(self.model.kmer_size, 1)

    def test_1_mer(self):
        all_1_mers = [''.join(x) for x in itertools.product('ACGT', repeat=1)]
        for kmer in all_1_mers:
            new_kmer = self.model.add_errors_to_kmer(kmer)
            self.assertEqual(new_kmer, [x for x in kmer])

    def test_2_mer(self):
        all_2_mers = [''.join(x) for x in itertools.product('ACGT', repeat=2)]
        for kmer in all_2_mers:
            new_kmer = self.model.add_errors_to_kmer(kmer)
            self.assertEqual(new_kmer, [x for x in kmer])

    def test_3_mer(self):
        all_3_mers = [''.join(x) for x in itertools.product('ACGT', repeat=3)]
        for kmer in all_3_mers:
            new_kmer = self.model.add_errors_to_kmer(kmer)
            self.assertEqual(new_kmer, [x for x in kmer])
