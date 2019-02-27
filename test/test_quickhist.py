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

import collections
import io
import unittest
import unittest.mock

import badread.quickhist


class TestHistogramFunctions(unittest.TestCase):

    def test_get_terminal_size(self):
        terminal_size = badread.quickhist.get_terminal_size_stderr()
        # The terminal could be any size, so just check that we're getting two integers.
        self.assertIsInstance(terminal_size[0], int)
        self.assertIsInstance(terminal_size[1], int)

    def test_get_max_width(self):
        self.assertTrue(80 <= badread.quickhist.get_max_width() <= 160)


FakeTerminalSize = collections.namedtuple('FakeTerminalSize', 'columns')


class TestGamma(unittest.TestCase):

    def setUp(self):
        self.captured_output = io.StringIO()
        self.a = 1.33136094675
        self.b = 0.0000887573964497
        self.n50 = 22200

    def tearDown(self):
        self.captured_output.close()

    def reset_output(self):
        self.tearDown()
        self.setUp()

    def test_quickhist_gamma(self):
        for height in [5, 10, 20]:
            badread.quickhist.quickhist_gamma(self.a, self.b, self.n50, height,
                                              output=self.captured_output)
            # The number of lines in the output should be twice the height (because there are two
            # histograms) plus 3 (2 for each x axis and 1 for the labels).
            self.assertEqual(self.captured_output.getvalue().count('\n'), 2 * height + 3)
            self.reset_output()

    def test_quickhist_gamma_small_terminal(self):
        with unittest.mock.patch('badread.quickhist.get_terminal_size_stderr',
                                 return_value=FakeTerminalSize(columns=100)):
            badread.quickhist.quickhist_gamma(self.a, self.b, self.n50, 10,
                                              output=self.captured_output)
        longest_line = max(len(x) for x in self.captured_output.getvalue().splitlines())
        self.assertLess(longest_line, 80)

    def test_quickhist_beta_gamma_terminal(self):
        with unittest.mock.patch('badread.quickhist.get_terminal_size_stderr',
                                 return_value=FakeTerminalSize(columns=140)):
            badread.quickhist.quickhist_gamma(self.a, self.b, self.n50, 10,
                                              output=self.captured_output)
        longest_line = max(len(x) for x in self.captured_output.getvalue().splitlines())
        self.assertGreater(longest_line, 80)


class TestBeta(unittest.TestCase):

    def setUp(self):
        self.captured_output = io.StringIO()
        self.a = 42.5
        self.b = 7.5
        self.max_identity = 100.0

    def tearDown(self):
        self.captured_output.close()

    def reset_output(self):
        self.tearDown()
        self.setUp()

    def test_quickhist_beta(self):
        for height in [5, 10, 20]:
            badread.quickhist.quickhist_beta(self.a, self.b, self.max_identity, height,
                                             output=self.captured_output)
            # The number of lines in the output should be the height plus 2 (1 for the x axis and 1
            # for the labels).
            self.assertEqual(self.captured_output.getvalue().count('\n'), height + 2)
            self.reset_output()

    def test_quickhist_beta_small_terminal(self):
        with unittest.mock.patch('badread.quickhist.get_terminal_size_stderr',
                                 return_value=FakeTerminalSize(columns=60)):
            badread.quickhist.quickhist_beta(self.a, self.b, self.max_identity, 10,
                                             output=self.captured_output)
        longest_line = max(len(x) for x in self.captured_output.getvalue().splitlines())
        self.assertLess(longest_line, 80)

    def test_quickhist_beta_big_terminal(self):
        with unittest.mock.patch('badread.quickhist.get_terminal_size_stderr',
                                 return_value=FakeTerminalSize(columns=180)):
            badread.quickhist.quickhist_beta(self.a, self.b, self.max_identity, 10,
                                             output=self.captured_output)
        longest_line = max(len(x) for x in self.captured_output.getvalue().splitlines())
        self.assertGreater(longest_line, 80)

    def test_quickhist_beta_exception(self):
        with unittest.mock.patch('os.get_terminal_size') as term_size_mock:
            term_size_mock.side_effect = AttributeError()
            badread.quickhist.quickhist_beta(self.a, self.b, self.max_identity, 10,
                                             output=self.captured_output)
