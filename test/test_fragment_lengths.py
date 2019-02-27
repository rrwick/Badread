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
import statistics
import unittest
import badread.fragment_lengths


class TestConstantFragmentLength(unittest.TestCase):

    def setUp(self):
        self.null = open(os.devnull, 'w')

    def tearDown(self):
        self.null.close()

    def test_constant_length_1(self):
        lengths = badread.fragment_lengths.FragmentLengths(1000, 0, output=self.null)
        for _ in range(100):
            self.assertEqual(lengths.get_fragment_length(), 1000)

    def test_constant_length_2(self):
        lengths = badread.fragment_lengths.FragmentLengths(5000, 0, output=self.null)
        for _ in range(100):
            self.assertEqual(lengths.get_fragment_length(), 5000)

    def test_constant_length_3(self):
        lengths = badread.fragment_lengths.FragmentLengths(8000, 0, output=self.null)
        for _ in range(100):
            self.assertEqual(lengths.get_fragment_length(), 8000)


class TestGammaFragmentLength(unittest.TestCase):

    def setUp(self):
        self.null = open(os.devnull, 'w')
        self.trials = 100000

    def tearDown(self):
        self.null.close()

    def test_gamma_length_1(self):
        lengths = badread.fragment_lengths.FragmentLengths(5000, 1000, output=self.null)
        all_lengths = [lengths.get_fragment_length() for _ in range(self.trials)]
        self.assertAlmostEqual(statistics.mean(all_lengths), 5000, delta=100)
        self.assertAlmostEqual(statistics.stdev(all_lengths), 1000, delta=100)

    def test_gamma_length_2(self):
        lengths = badread.fragment_lengths.FragmentLengths(5000, 3000, output=self.null)
        all_lengths = [lengths.get_fragment_length() for _ in range(self.trials)]
        self.assertAlmostEqual(statistics.mean(all_lengths), 5000, delta=100)
        self.assertAlmostEqual(statistics.stdev(all_lengths), 3000, delta=100)

    def test_gamma_length_3(self):
        lengths = badread.fragment_lengths.FragmentLengths(20000, 30000, output=self.null)
        all_lengths = [lengths.get_fragment_length() for _ in range(self.trials)]
        self.assertAlmostEqual(statistics.mean(all_lengths), 20000, delta=1000)
        self.assertAlmostEqual(statistics.stdev(all_lengths), 30000, delta=1000)
