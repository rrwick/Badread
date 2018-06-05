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

import edlib
import os
import pathlib
import random
import unittest
import badread.simulate
import badread.identities
import badread.error_model
import badread.misc


class TestPerfectSequenceFragment(unittest.TestCase):
    """
    Tests the sequence_fragment function with a perfect error model.
    """
    def setUp(self):
        self.null = open(os.devnull, 'w')
        self.model = badread.error_model.ErrorModel('perfect', output=self.null)

    def tearDown(self):
        self.null.close()

    def test_sequence_fragment_1(self):
        frag = 'GACCCAGTTTTTTTACTGATTCAGCGTAGGTGCTCTGATCTTCACGCATCTTTGACCGCC'
        seq, qual = badread.simulate.sequence_fragment(frag, 1.0, None, None, self.model)
        self.assertEqual(frag, seq)
        self.assertEqual(len(frag), len(qual))
        for q in qual:
            self.assertTrue(q in 'ABCDEFGHI')

    def test_sequence_fragment_2(self):
        # Beta distribution parameters are ignored for perfect error models.
        frag = 'TATAAAGACCCCACTTTTGAAGCCAGAGGTAATGGCCGTGATGGCGTTAAATTCCCTTCC'
        seq, qual = badread.simulate.sequence_fragment(frag, 0.9, None, None, self.model)
        self.assertEqual(frag, seq)
        self.assertEqual(len(frag), len(qual))
        for q in qual:
            self.assertTrue(q in 'ABCDEFGHI')


class TestRandomErrorSequenceFragment(unittest.TestCase):
    """
    Tests the sequence_fragment function with a random error model.
    """
    def setUp(self):
        self.null = open(os.devnull, 'w')
        self.model = badread.error_model.ErrorModel('random', output=self.null)
        self.trials = 100

    def tearDown(self):
        self.null.close()

    def test_random_70_identity(self):
        identity_sum = 0.0
        for _ in range(self.trials):
            frag = badread.misc.get_random_sequence(random.randint(50, 1000))
            seq, qual = badread.simulate.sequence_fragment(frag, 0.7, None, None, self.model)
            identity_sum += 1.0 - (edlib.align(frag, seq, task='path')['editDistance'] / len(frag))
        mean_identity = identity_sum / self.trials
        self.assertAlmostEqual(mean_identity, 0.7, delta=0.001)

    def test_random_80_identity(self):
        identity_sum = 0.0
        for _ in range(self.trials):
            frag = badread.misc.get_random_sequence(random.randint(50, 1000))
            seq, qual = badread.simulate.sequence_fragment(frag, 0.8, None, None, self.model)
            cigar = edlib.align(frag, seq, task='path')['cigar']
            identity_sum += badread.error_model.identity_from_edlib_cigar(cigar)
        mean_identity = identity_sum / self.trials
        self.assertAlmostEqual(mean_identity, 0.8, delta=0.001)

    def test_random_90_identity(self):
        identity_sum = 0.0
        for _ in range(self.trials):
            frag = badread.misc.get_random_sequence(random.randint(50, 1000))
            seq, qual = badread.simulate.sequence_fragment(frag, 0.9, None, None, self.model)
            cigar = edlib.align(frag, seq, task='path')['cigar']
            identity_sum += badread.error_model.identity_from_edlib_cigar(cigar)
        mean_identity = identity_sum / self.trials
        self.assertAlmostEqual(mean_identity, 0.9, delta=0.001)

    def test_random_100_identity(self):
        identity_sum = 0.0
        for _ in range(self.trials):
            frag = badread.misc.get_random_sequence(random.randint(50, 1000))
            seq, qual = badread.simulate.sequence_fragment(frag, 1.0, None, None, self.model)
            cigar = edlib.align(frag, seq, task='path')['cigar']
            identity_sum += badread.error_model.identity_from_edlib_cigar(cigar)
        mean_identity = identity_sum / self.trials
        self.assertEqual(mean_identity, 1.0)


class TestNanoporeSequenceFragment(unittest.TestCase):
    """
    Tests the sequence_fragment function with a Nanopore-based error model.
    """
    def setUp(self):
        self.null = open(os.devnull, 'w')
        self.trials = 100
        model_file = pathlib.Path(__file__).parent.parent / 'error_models' / 'nanopore_7-mer_model'
        self.model = badread.error_model.ErrorModel(model_file, output=self.null)

    def tearDown(self):
        self.null.close()
    #
    # def test_nanopore_70_identity(self):
    #     identity_sum = 0.0
    #     for _ in range(self.trials):
    #         frag = badread.misc.get_random_sequence(random.randint(50, 1000))
    #         seq, qual = badread.simulate.sequence_fragment(frag, 0.7, None, None, self.model)
    #         cigar = edlib.align(frag, seq, task='path')['cigar']
    #         identity_sum += badread.error_model.identity_from_edlib_cigar(cigar)
    #     mean_identity = identity_sum / self.trials
    #     self.assertAlmostEqual(mean_identity, 0.7, delta=0.001)
    #
    # def test_nanopore_80_identity(self):
    #     identity_sum = 0.0
    #     for _ in range(self.trials):
    #         frag = badread.misc.get_random_sequence(random.randint(50, 1000))
    #         seq, qual = badread.simulate.sequence_fragment(frag, 0.8, None, None, self.model)
    #         cigar = edlib.align(frag, seq, task='path')['cigar']
    #         identity_sum += badread.error_model.identity_from_edlib_cigar(cigar)
    #     mean_identity = identity_sum / self.trials
    #     self.assertAlmostEqual(mean_identity, 0.8, delta=0.001)

    def test_nanopore_90_identity(self):
        identity_sum = 0.0
        for _ in range(self.trials):
            frag = badread.misc.get_random_sequence(random.randint(50, 1000))
            seq, qual = badread.simulate.sequence_fragment(frag, 0.9, None, None, self.model)
            cigar = edlib.align(frag, seq, task='path')['cigar']
            identity_sum += badread.error_model.identity_from_edlib_cigar(cigar)
        mean_identity = identity_sum / self.trials
        self.assertAlmostEqual(mean_identity, 0.9, delta=0.001)

    # def test_nanopore_100_identity(self):
    #     identity_sum = 0.0
    #     for _ in range(self.trials):
    #         frag = badread.misc.get_random_sequence(random.randint(50, 1000))
    #         seq, qual = badread.simulate.sequence_fragment(frag, 1.0, None, None, self.model)
    #         cigar = edlib.align(frag, seq, task='path')['cigar']
    #         identity_sum += badread.error_model.identity_from_edlib_cigar(cigar)
    #     mean_identity = identity_sum / self.trials
    #     self.assertEqual(mean_identity, 1.0)
