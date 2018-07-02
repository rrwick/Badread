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

import unittest

import badread.fragment_lengths
import badread.misc
import badread.simulate


class TestParameterParsing(unittest.TestCase):

    def test_params_1(self):
        rate, amount = badread.simulate.adapter_parameters('10,60')
        self.assertEqual(rate, 0.1)
        self.assertEqual(amount, 0.6)

    def test_params_2(self):
        rate, amount = badread.simulate.adapter_parameters('0,0')
        self.assertEqual(rate, 0.0)
        self.assertEqual(amount, 0.0)

    def test_params_3(self):
        with self.assertRaises(SystemExit):
            _, _ = badread.simulate.adapter_parameters('1')

    def test_params_4(self):
        with self.assertRaises(SystemExit):
            _, _ = badread.simulate.adapter_parameters('wrong_input')

    def test_params_5(self):
        with self.assertRaises(SystemExit):
            _, _ = badread.simulate.adapter_parameters('6,bad')

    def test_params_6(self):
        with self.assertRaises(SystemExit):
            _, _ = badread.simulate.adapter_parameters('bad,10')


class TestStartAdapters(unittest.TestCase):

    def setUp(self):
        self.trials = 100
        self.seq = 'TGTATAATACGACGCCGAGC'

    def test_start_adapters_1(self):
        # A rate and amount of 1 means the adapter is always returned in its entirety.
        rate, amount = 1.0, 1.0
        for _ in range(self.trials):
            adapter = badread.simulate.get_start_adapter(rate, amount, self.seq)
            self.assertEqual(adapter, self.seq)

    def test_start_adapters_2(self):
        # A rate of 0.5 and amount of 1 means we either get nothing or the full adapter.
        rate, amount = 0.5, 1.0
        results = set()
        for _ in range(self.trials):
            results.add(badread.simulate.get_start_adapter(rate, amount, self.seq))
        self.assertEqual(len(results), 2)
        self.assertTrue('' in results)
        self.assertTrue(self.seq in results)

    def test_start_adapters_3(self):
        # A rate of 1.0 and amount of 0.5 means we will get all sorts of lengths of adapter.
        rate, amount = 1.0, 0.5
        lengths = set()
        for _ in range(self.trials):
            adapter = badread.simulate.get_start_adapter(rate, amount, self.seq)
            lengths.add(len(adapter))
            self.assertTrue(self.seq.endswith(adapter))
        self.assertGreater(len(lengths), 15)
        self.assertLessEqual(max(lengths), len(self.seq))

    def test_no_start_adapters_1(self):
        rate, amount = 0.0, 1.0
        for _ in range(self.trials):
            adapter = badread.simulate.get_start_adapter(rate, amount, self.seq)
            self.assertEqual(adapter, '')

    def test_no_start_adapters_2(self):
        rate, amount = 1.0, 0.0
        for _ in range(self.trials):
            adapter = badread.simulate.get_start_adapter(rate, amount, self.seq)
            self.assertEqual(adapter, '')

    def test_no_start_adapters_3(self):
        rate, amount = 1.0, 1.0
        for _ in range(self.trials):
            adapter = badread.simulate.get_start_adapter(rate, amount, '')
            self.assertEqual(adapter, '')


class TestEndAdapters(unittest.TestCase):

    def setUp(self):
        self.trials = 100
        self.seq = 'ATAACAAACGCTAATTGCAA'

    def test_start_adapters_1(self):
        # A rate and amount of 1 means the adapter is always returned in its entirety.
        rate, amount = 1.0, 1.0
        for _ in range(self.trials):
            adapter = badread.simulate.get_end_adapter(rate, amount, self.seq)
            self.assertEqual(adapter, self.seq)

    def test_start_adapters_2(self):
        # A rate of 0.5 and amount of 1 means we either get nothing or the full adapter.
        rate, amount = 0.5, 1.0
        results = set()
        for _ in range(self.trials):
            results.add(badread.simulate.get_end_adapter(rate, amount, self.seq))
        self.assertEqual(len(results), 2)
        self.assertTrue('' in results)
        self.assertTrue(self.seq in results)

    def test_end_adapters_3(self):
        # A rate of 1.0 and amount of 0.5 means we will get all sorts of lengths of adapter.
        rate, amount = 1.0, 0.5
        lengths = set()
        for _ in range(self.trials):
            adapter = badread.simulate.get_end_adapter(rate, amount, self.seq)
            lengths.add(len(adapter))
            self.assertTrue(self.seq.startswith(adapter))
        self.assertGreater(len(lengths), 15)
        self.assertLessEqual(max(lengths), len(self.seq))

    def test_no_end_adapters_1(self):
        rate, amount = 0.0, 1.0
        for _ in range(self.trials):
            adapter = badread.simulate.get_end_adapter(rate, amount, self.seq)
            self.assertEqual(adapter, '')

    def test_no_end_adapters_2(self):
        rate, amount = 1.0, 0.0
        for _ in range(self.trials):
            adapter = badread.simulate.get_end_adapter(rate, amount, self.seq)
            self.assertEqual(adapter, '')

    def test_no_end_adapters_3(self):
        rate, amount = 1.0, 1.0
        for _ in range(self.trials):
            adapter = badread.simulate.get_end_adapter(rate, amount, '')
            self.assertEqual(adapter, '')
