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

import os
import unittest
import badread.error_model
import badread.identities


class TestConstantIdentity(unittest.TestCase):

    def setUp(self):
        self.null = open(os.devnull, 'w')

    def tearDown(self):
        self.null.close()

    def test_constant_identity_1(self):
        identities = badread.identities.Identities('constant', 100, output=self.null)
        for _ in range(100):
            self.assertEqual(identities.get_identity(), 1.0)

    def test_constant_identity_2(self):
        identities = badread.identities.Identities('constant', 80, output=self.null)
        for _ in range(100):
            self.assertEqual(identities.get_identity(), 0.8)


class TestBetaIdentity(unittest.TestCase):

    def setUp(self):
        self.null = open(os.devnull, 'w')
        self.trials = 100000

    def tearDown(self):
        self.null.close()

    def test_beta_identity_1(self):
        identities = badread.identities.Identities('beta', 90, 4, 100, output=self.null)
        mean = sum(identities.get_identity() for _ in range(self.trials)) / self.trials
        self.assertAlmostEqual(mean, 0.9, delta=0.01)

    def test_beta_identity_2(self):
        identities = badread.identities.Identities('beta', 90, 4, 95, output=self.null)
        mean = sum(identities.get_identity() for _ in range(self.trials)) / self.trials
        self.assertAlmostEqual(mean, 0.9, delta=0.01)

    def test_beta_identity_3(self):
        identities = badread.identities.Identities('beta', 90, 4, 90, output=self.null)
        mean = sum(identities.get_identity() for _ in range(self.trials)) / self.trials
        self.assertAlmostEqual(mean, 0.9, delta=0.01)

    def test_beta_identity_4(self):
        identities = badread.identities.Identities('beta', 90, 3, 100, output=self.null)
        mean = sum(identities.get_identity() for _ in range(self.trials)) / self.trials
        self.assertAlmostEqual(mean, 0.9, delta=0.01)

    def test_beta_identity_5(self):
        identities = badread.identities.Identities('beta', 90, 2, 100, output=self.null)
        mean = sum(identities.get_identity() for _ in range(self.trials)) / self.trials
        self.assertAlmostEqual(mean, 0.9, delta=0.01)

    def test_beta_identity_6(self):
        identities = badread.identities.Identities('beta', 90, 8, 100, output=self.null)
        mean = sum(identities.get_identity() for _ in range(self.trials)) / self.trials
        self.assertAlmostEqual(mean, 0.9, delta=0.01)


class TestPerfectErrorModel(unittest.TestCase):
    """
    If the error model was set to 'perfect', then we always use 100% identity, regardless of the
    settings.
    """

    def setUp(self):
        self.null = open(os.devnull, 'w')
        self.model = badread.error_model.ErrorModel('perfect', output=self.null)

    def tearDown(self):
        self.null.close()

    def test_perfect_error_model_1(self):
        identities = badread.identities.Identities('constant', 100, error_model=self.model,
                                                   output=self.null)
        for _ in range(10):
            self.assertEqual(identities.get_identity(), 1.0)

    def test_perfect_error_model_2(self):
        identities = badread.identities.Identities('constant', 80, error_model=self.model,
                                                   output=self.null)
        for _ in range(10):
            self.assertEqual(identities.get_identity(), 1.0)

    def test_perfect_error_model_3(self):
        identities = badread.identities.Identities('beta', 90, 4, 95, error_model=self.model,
                                                   output=self.null)
        for _ in range(10):
            self.assertEqual(identities.get_identity(), 1.0)
