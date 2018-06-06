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


class TestSequenceFragment(unittest.TestCase):
    """
    Tests the sequence_fragment function with a random error model.
    """
    def setUp(self):
        self.null = open(os.devnull, 'w')
        self.trials = 20
        self.identities_to_test = [1.0, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65]
        self.read_lengths_to_test = [10000, 1000]
        self.read_delta = 0.01
        self.mean_delta = 0.002
        print()  # TEMP

    def tearDown(self):
        self.null.close()

    def identity_test(self, target_identity, read_length, model):
        identity_sum = 0.0
        for _ in range(self.trials):
            frag = badread.misc.get_random_sequence(read_length)
            seq, qual = badread.simulate.sequence_fragment(frag, target_identity, None, None, model)
            cigar = edlib.align(frag, seq, task='path')['cigar']
            read_identity = badread.error_model.identity_from_edlib_cigar(cigar)
            identity_sum += read_identity

            # if target_identity == 1.0:
            #     self.assertEqual(read_identity, 1.0)
            # else:
            #     self.assertAlmostEqual(read_identity, target_identity, delta=self.read_delta)

        mean_identity = identity_sum / self.trials
        # if target_identity == 1.0:
        #     self.assertEqual(mean_identity, 1.0)
        # else:
        #     self.assertAlmostEqual(mean_identity, target_identity, delta=self.mean_delta)

        print('\t'.join(str(x) for x in [read_length, model.type, target_identity, mean_identity]))  # TEMP

    def test_random_identity(self):
        model = badread.error_model.ErrorModel('random', output=self.null)
        for identity in self.identities_to_test:
            for read_length in self.read_lengths_to_test:
                self.identity_test(identity, read_length, model)

    def test_nanopore_identity(self):
        model_file = pathlib.Path(__file__).parent.parent / 'error_models' / 'nanopore_7-mer_model'
        model = badread.error_model.ErrorModel(model_file, output=self.null)
        for identity in self.identities_to_test:
            for read_length in self.read_lengths_to_test:
                self.identity_test(identity, read_length, model)

    def test_pacbio_identity(self):
        model_file = pathlib.Path(__file__).parent.parent / 'error_models' / 'pacbio_7-mer_model'
        model = badread.error_model.ErrorModel(model_file, output=self.null)
        for identity in self.identities_to_test:
            for read_length in self.read_lengths_to_test:
                self.identity_test(identity, read_length, model)
