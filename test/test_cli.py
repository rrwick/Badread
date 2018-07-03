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

from contextlib import contextmanager
from io import StringIO
import os
import sys
import unittest

import badread.badread


@contextmanager
def captured_output():
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err


class TestCommandLine(unittest.TestCase):

    def setUp(self):
        self.ref = os.path.join(os.path.dirname(__file__), 'test_ref_1.fasta')

    def check_simulate_error(self, args, message, ref=None):
        if ref is None:
            ref = self.ref
        args = ['simulate', '--reference', ref, '--quantity', '50x'] + args
        args = badread.badread.parse_args(args)
        with self.assertRaises(SystemExit) as cm:
            badread.badread.check_simulate_args(args)
        self.assertNotEqual(cm.exception.code, 0)
        self.assertTrue(message in str(cm.exception))

    def test_no_args(self):
        # When called with no arguments, the program behaves like it's hit an error: the help text
        # goes to stderr and it returns a non-zero exit code.
        with captured_output() as (out, err):
            with self.assertRaises(SystemExit) as cm:
                badread.badread.parse_args([])
        self.assertNotEqual(cm.exception.code, 0)
        self.assertTrue('Badread' in err.getvalue().strip())

    def test_help_1(self):
        # When called with -h or --help, the program doesn't behaves like it's hit an error: the
        # help text goes to stdout and it returns a exit code of zero.
        with captured_output() as (out, err):
            with self.assertRaises(SystemExit) as cm:
                badread.badread.parse_args(['-h'])
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(err.getvalue().strip(), '')
        self.assertTrue('Badread' in out.getvalue().strip())

    def test_help_2(self):
        # When called with -h or --help, the program doesn't behave like it's hit an error: the
        # help text goes to stdout and it returns a exit code of zero.
        with captured_output() as (out, err):
            with self.assertRaises(SystemExit) as cm:
                badread.badread.parse_args(['--help'])
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(err.getvalue().strip(), '')
        self.assertTrue('Badread' in out.getvalue().strip())

    def test_simulate_args(self):
        args = badread.badread.parse_args(['simulate', '--reference', self.ref,
                                           '--quantity', '50x', '--length', '5432,123',
                                           '--identity', '78,89,4', '--glitches', '9876,12,34'])
        self.assertEqual(args.reference, self.ref)
        self.assertEqual(args.quantity, '50x')
        badread.badread.check_simulate_args(args)
        self.assertEqual(args.mean_frag_length, 5432)
        self.assertEqual(args.frag_length_stdev, 123)
        self.assertEqual(args.mean_identity, 78)
        self.assertEqual(args.max_identity, 89)
        self.assertEqual(args.identity_stdev, 4)
        self.assertEqual(args.glitch_rate, 9876)
        self.assertEqual(args.glitch_size, 12)
        self.assertEqual(args.glitch_skip, 34)

    def test_simulate_help(self):
        with captured_output() as (out, err):
            with self.assertRaises(SystemExit) as cm:
                badread.badread.parse_args(['simulate'])
        self.assertNotEqual(cm.exception.code, 0)
        self.assertEqual(out.getvalue().strip(), '')
        self.assertTrue('simulate' in err.getvalue().strip())

    def test_bad_simulated_args_1(self):
        self.check_simulate_error(['--chimera', '101'], 'cannot be greater than')

    def test_bad_simulated_args_2(self):
        self.check_simulate_error(['--junk_reads', '101'], 'cannot be greater than')

    def test_bad_simulated_args_3(self):
        self.check_simulate_error(['--random_reads', '101'], 'cannot be greater than')

    def test_bad_simulated_args_4(self):
        self.check_simulate_error(['--junk_reads', '60', '--random_reads', '60'],
                                  'cannot sum to more than')

    def test_bad_simulated_args_5(self):
        self.check_simulate_error(['--length', '1,1'], 'must be at least')

    def test_bad_simulated_args_6(self):
        self.check_simulate_error(['--length', '1000,-1'], 'cannot be negative')

    def test_bad_simulated_args_7(self):
        self.check_simulate_error(['--length', 'abc'], 'could not parse')

    def test_bad_simulated_args_8(self):
        self.check_simulate_error(['--length', '6,A'], 'could not parse')

    def test_bad_simulated_args_9(self):
        self.check_simulate_error(['--identity', 'abc'], 'could not parse')

    def test_bad_simulated_args_10(self):
        self.check_simulate_error(['--identity', '85,A,5'], 'could not parse')

    def test_bad_simulated_args_11(self):
        self.check_simulate_error(['--identity', '85,101,5'], 'cannot be more than')

    def test_bad_simulated_args_12(self):
        self.check_simulate_error(['--identity', '101,102,5'], 'cannot be more than')

    def test_bad_simulated_args_13(self):
        self.check_simulate_error(['--identity', '85,80,5'], 'cannot be larger than max')

    def test_bad_simulated_args_14(self):
        self.check_simulate_error(['--identity', '85,90,-3'], 'cannot be negative')

    def test_bad_simulated_args_15(self):
        self.check_simulate_error(['--identity', '30,90,5'], 'must be at least')

    def test_bad_simulated_args_16(self):
        self.check_simulate_error(['--identity', '90,30,5'], 'must be at least')

    def test_bad_simulated_args_17(self):
        self.check_simulate_error(['--glitches', 'abc'], 'could not parse')

    def test_bad_simulated_args_18(self):
        self.check_simulate_error(['--glitches', '500,Z,10'], 'could not parse')

    def test_bad_simulated_args_19(self):
        self.check_simulate_error(['--glitches', '500,-10,10'], 'must contain non-negative')

    def test_bad_simulated_args_20(self):
        self.check_simulate_error(['--error_model', 'not_a_file'], 'or a filename')

    def test_bad_simulated_args_21(self):
        self.check_simulate_error(['--qscore_model', 'not_a_file'], 'or a filename')

    def test_bad_simulated_args_22(self):
        self.check_simulate_error([], 'not a file', ref='not_a_file')
