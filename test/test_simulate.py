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

import edlib
import os
import pathlib
import statistics
import unittest

import badread.simulate
import badread.identities
import badread.error_model
import badread.qscore_model
import badread.misc


VERBOSE = False  # Turn this on to see detailed read identity output


class TestPerfectSequenceFragment(unittest.TestCase):
    """
    Tests the sequence_fragment function with a perfect error model.
    """
    def setUp(self):
        self.null = open(os.devnull, 'w')
        self.error_model = badread.error_model.ErrorModel('random', output=self.null)
        self.qscore_model = badread.qscore_model.QScoreModel('random', output=self.null)

    def tearDown(self):
        self.null.close()

    def test_perfect_sequence_fragment(self):
        frag = 'GACCCAGTTTTTTTACTGATTCAGCGTAGGTGCTCTGATCTTCACGCATCTTTGACCGCC'
        seq, qual, _, _ = badread.simulate.sequence_fragment(frag, 1.0, self.error_model,
                                                             self.qscore_model)
        self.assertEqual(frag, seq)
        self.assertEqual(len(frag), len(qual))


class TestSequenceFragment(unittest.TestCase):
    """
    Tests the sequence_fragment function with a random error model.
    """
    def setUp(self):
        self.null = open(os.devnull, 'w')
        self.trials = 20
        self.identities_to_test = [1.0, 0.9, 0.8]
        self.read_lengths_to_test = [3000, 1000]
        self.read_delta = 0.5
        self.mean_delta = 0.05
        self.repo_dir = pathlib.Path(__file__).parent.parent

    def tearDown(self):
        self.null.close()

    def identity_test(self, target_identity, read_length, error_model, qscore_model):
        target_errors = 1.0 - target_identity
        read_delta = self.read_delta * target_errors
        mean_delta = self.mean_delta * target_errors

        if VERBOSE:
            print(f'\nRead length: {read_length}, target identity: {target_identity}')
            print(f'    allowed error per read:   {read_delta:.4f}')
            print(f'    allowed error in mean:    {mean_delta:.4f}')
            print('    identities: ', end='')

        read_identities = []
        for i in range(self.trials):
            frag = badread.misc.get_random_sequence(read_length)
            seq, qual, _, _ = badread.simulate.sequence_fragment(frag, target_identity,
                                                                 error_model, qscore_model)
            cigar = edlib.align(frag, seq, task='path')['cigar']
            read_identity = badread.misc.identity_from_edlib_cigar(cigar)
            read_identities.append(read_identity)

            if VERBOSE:
                print('{:.4f}'.format(read_identity), flush=True,
                      end='\n                ' if (i+1) % 20 == 0 else ' ')

            self.assertAlmostEqual(read_identity, target_identity, delta=read_delta)

        mean_identity = statistics.mean(read_identities)
        if VERBOSE:
            print('\r' if self.trials % 20 == 0 else '\n', end='')
            print(f'    mean:       {mean_identity:.4f}')

        self.assertAlmostEqual(mean_identity, target_identity, delta=mean_delta)
        if VERBOSE:
            print('    PASS')

    def test_random_identity(self):
        if VERBOSE:
            print('\n\nRANDOM ERROR MODEL\n------------------')
        error_model = badread.error_model.ErrorModel('random', output=self.null)
        qscore_model = badread.qscore_model.QScoreModel('random', output=self.null)
        for identity in self.identities_to_test:
            for read_length in self.read_lengths_to_test:
                self.identity_test(identity, read_length, error_model, qscore_model)

    def test_nanopore2018_identity(self):
        if VERBOSE:
            print('\n\nNANOPORE ERROR MODEL\n--------------------')
        model_file = self.repo_dir / 'badread' / 'error_models' / 'nanopore2018.gz'
        error_model = badread.error_model.ErrorModel(model_file, output=self.null)
        qscore_model = badread.qscore_model.QScoreModel('random', output=self.null)
        for identity in self.identities_to_test:
            for read_length in self.read_lengths_to_test:
                self.identity_test(identity, read_length, error_model, qscore_model)

    def test_nanopore2020_identity(self):
        if VERBOSE:
            print('\n\nNANOPORE ERROR MODEL\n--------------------')
        model_file = self.repo_dir / 'badread' / 'error_models' / 'nanopore2020.gz'
        error_model = badread.error_model.ErrorModel(model_file, output=self.null)
        qscore_model = badread.qscore_model.QScoreModel('random', output=self.null)
        for identity in self.identities_to_test:
            for read_length in self.read_lengths_to_test:
                self.identity_test(identity, read_length, error_model, qscore_model)

    def test_nanopore2023_identity(self):
        if VERBOSE:
            print('\n\nNANOPORE ERROR MODEL\n--------------------')
        model_file = self.repo_dir / 'badread' / 'error_models' / 'nanopore2023.gz'
        error_model = badread.error_model.ErrorModel(model_file, output=self.null)
        qscore_model = badread.qscore_model.QScoreModel('random', output=self.null)
        for identity in self.identities_to_test:
            for read_length in self.read_lengths_to_test:
                self.identity_test(identity, read_length, error_model, qscore_model)

    def test_pacbio2016_identity(self):
        if VERBOSE:
            print('\n\nPACBIO ERROR MODEL\n------------------')
        model_file = self.repo_dir / 'badread' / 'error_models' / 'pacbio2016.gz'
        error_model = badread.error_model.ErrorModel(model_file, output=self.null)
        qscore_model = badread.qscore_model.QScoreModel('random', output=self.null)
        for identity in self.identities_to_test:
            for read_length in self.read_lengths_to_test:
                self.identity_test(identity, read_length, error_model, qscore_model)

    def test_pacbio2021_identity(self):
        if VERBOSE:
            print('\n\nPACBIO ERROR MODEL\n------------------')
        model_file = self.repo_dir / 'badread' / 'error_models' / 'pacbio2021.gz'
        error_model = badread.error_model.ErrorModel(model_file, output=self.null)
        qscore_model = badread.qscore_model.QScoreModel('random', output=self.null)
        for identity in self.identities_to_test:
            for read_length in self.read_lengths_to_test:
                self.identity_test(identity, read_length, error_model, qscore_model)
