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
import badread.error_model
import badread.identities


class TestConstantIdentity(unittest.TestCase):

    def setUp(self):
        self.null = open(os.devnull, 'w')
        self.trials = 20

    def tearDown(self):
        self.null.close()

    def test_constant_identity_1(self):
        identities = badread.identities.Identities(100, 4, 100, output=self.null)
        for _ in range(self.trials):
            self.assertEqual(identities.get_identity(), 1.0)

    def test_constant_identity_2(self):
        identities = badread.identities.Identities(80, 4, 80, output=self.null)
        for _ in range(self.trials):
            self.assertEqual(identities.get_identity(), 0.8)

    def test_constant_identity_3(self):
        identities = badread.identities.Identities(90, 0, 100, output=self.null)
        for _ in range(self.trials):
            self.assertEqual(identities.get_identity(), 0.9)


class TestBetaIdentity(unittest.TestCase):

    def setUp(self):
        self.null = open(os.devnull, 'w')
        self.trials = 100000

    def tearDown(self):
        self.null.close()

    def test_beta_identity_1(self):
        identities = badread.identities.Identities(90, 4, 100, output=self.null)
        mean = sum(identities.get_identity() for _ in range(self.trials)) / self.trials
        self.assertAlmostEqual(mean, 0.9, delta=0.01)

    def test_beta_identity_2(self):
        identities = badread.identities.Identities(90, 4, 95, output=self.null)
        mean = sum(identities.get_identity() for _ in range(self.trials)) / self.trials
        self.assertAlmostEqual(mean, 0.9, delta=0.01)

    def test_beta_identity_3(self):
        identities = badread.identities.Identities(90, 4, 90, output=self.null)
        mean = sum(identities.get_identity() for _ in range(self.trials)) / self.trials
        self.assertAlmostEqual(mean, 0.9, delta=0.01)

    def test_beta_identity_4(self):
        identities = badread.identities.Identities(90, 3, 100, output=self.null)
        mean = sum(identities.get_identity() for _ in range(self.trials)) / self.trials
        self.assertAlmostEqual(mean, 0.9, delta=0.01)

    def test_beta_identity_5(self):
        identities = badread.identities.Identities(90, 2, 100, output=self.null)
        mean = sum(identities.get_identity() for _ in range(self.trials)) / self.trials
        self.assertAlmostEqual(mean, 0.9, delta=0.01)

    def test_beta_identity_6(self):
        identities = badread.identities.Identities(90, 8, 100, output=self.null)
        mean = sum(identities.get_identity() for _ in range(self.trials)) / self.trials
        self.assertAlmostEqual(mean, 0.9, delta=0.01)

    def test_bad_identity(self):
        with self.assertRaises(SystemExit) as cm:
            identities = badread.identities.Identities(81.9, 5.5, 82.1, output=self.null)
            identities.get_identity()
        self.assertTrue('invalid beta parameters' in str(cm.exception))


class TestNormalIdentity(unittest.TestCase):

    def setUp(self):
        self.null = open(os.devnull, 'w')
        self.trials = 100000

    def tearDown(self):
        self.null.close()

    def test_normal_identity_1(self):
        identities = badread.identities.Identities(20, 2, None, output=self.null)
        mean = sum(identities.get_identity() for _ in range(self.trials)) / self.trials
        self.assertAlmostEqual(mean, 0.98888, delta=0.01)

    def test_normal_identity_2(self):
        identities = badread.identities.Identities(10, 0, None, output=self.null)
        mean = sum(identities.get_identity() for _ in range(self.trials)) / self.trials
        self.assertAlmostEqual(mean, 0.9, delta=0.01)

    def test_normal_identity_3(self):
        identities = badread.identities.Identities(20, 0, None, output=self.null)
        mean = sum(identities.get_identity() for _ in range(self.trials)) / self.trials
        self.assertAlmostEqual(mean, 0.99, delta=0.01)

    def test_normal_identity_4(self):
        identities = badread.identities.Identities(30, 0, None, output=self.null)
        mean = sum(identities.get_identity() for _ in range(self.trials)) / self.trials
        self.assertAlmostEqual(mean, 0.999, delta=0.01)
