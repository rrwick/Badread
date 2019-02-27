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

import unittest

import badread.simulate


class TestTargetSize(unittest.TestCase):

    def test_relative(self):
        self.assertEqual(badread.simulate.get_target_size(1000000, '0x'), 0)
        self.assertEqual(badread.simulate.get_target_size(1000000, '0.0x'), 0)
        self.assertEqual(badread.simulate.get_target_size(1000000, '0.5x'), 500000)
        self.assertEqual(badread.simulate.get_target_size(1000000, '0.50000x'), 500000)
        self.assertEqual(badread.simulate.get_target_size(1000000, '0000.5x'), 500000)
        self.assertEqual(badread.simulate.get_target_size(1000000, '1x'), 1000000)
        self.assertEqual(badread.simulate.get_target_size(1000000, '1.0x'), 1000000)
        self.assertEqual(badread.simulate.get_target_size(1000000, '2x'), 2000000)
        self.assertEqual(badread.simulate.get_target_size(1000000, '2.5x'), 2500000)
        self.assertEqual(badread.simulate.get_target_size(1000000, '10x'), 10000000)
        self.assertEqual(badread.simulate.get_target_size(1000000, '100x'), 100000000)
        self.assertEqual(badread.simulate.get_target_size(1000000, '1000x'), 1000000000)
        self.assertEqual(badread.simulate.get_target_size(1000000, '1000.1x'), 1000100000)

    def test_absolute(self):
        self.assertEqual(badread.simulate.get_target_size(1000000, '0'), 0)
        self.assertEqual(badread.simulate.get_target_size(1000000, '10'), 10)
        self.assertEqual(badread.simulate.get_target_size(1000000, '237864273'), 237864273)

    def test_kilo(self):
        self.assertEqual(badread.simulate.get_target_size(1000000, '0k'), 0)
        self.assertEqual(badread.simulate.get_target_size(1000000, '10k'), 10000)
        self.assertEqual(badread.simulate.get_target_size(1000000, '5678K'), 5678000)

    def test_mega(self):
        self.assertEqual(badread.simulate.get_target_size(1000000, '0M'), 0)
        self.assertEqual(badread.simulate.get_target_size(1000000, '10.01m'), 10010000)
        self.assertEqual(badread.simulate.get_target_size(1000000, '12345m'), 12345000000)

    def test_giga(self):
        self.assertEqual(badread.simulate.get_target_size(1000000, '0g'), 0)
        self.assertEqual(badread.simulate.get_target_size(1000000, '3.1g'), 3100000000)
        self.assertEqual(badread.simulate.get_target_size(1000000, '987654321G'),
                         987654321000000000)

    def test_round(self):
        self.assertEqual(badread.simulate.get_target_size(1000000, '10.0001k'), 10000)
        self.assertEqual(badread.simulate.get_target_size(1000000, '10.0009k'), 10001)
        self.assertEqual(badread.simulate.get_target_size(10, '0.1x'), 1)
        self.assertEqual(badread.simulate.get_target_size(10, '0.14x'), 1)
        self.assertEqual(badread.simulate.get_target_size(10, '0.16x'), 2)

    def test_bad_format(self):
        with self.assertRaises(SystemExit):
            badread.simulate.get_target_size(1000000, '0T')
        with self.assertRaises(SystemExit):
            badread.simulate.get_target_size(1000000, 'abcdefg')
        with self.assertRaises(SystemExit):
            badread.simulate.get_target_size(1000000, '10.0a0k')
