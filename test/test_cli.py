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
import unittest.mock
import sys

import badread.__main__
import badread.misc


class TestWholeCommands(unittest.TestCase):

    def setUp(self):
        self.ref_filename = os.path.join(os.path.dirname(__file__), 'test_alignment_ref.fasta')
        self.reads_filename = os.path.join(os.path.dirname(__file__), 'test_alignment_reads.fastq')
        self.paf_filename = os.path.join(os.path.dirname(__file__), 'test_alignment.paf')
        self.null = open(os.devnull, 'w')

    def tearDown(self):
        self.null.close()

    def test_simulate(self):
        test_args = ['badread', 'simulate', '--reference', self.ref_filename, '--quantity', '1x',
                     '--error_model', 'random', '--qscore_model', 'random']
        with unittest.mock.patch.object(sys, 'argv', test_args):
            with badread.misc.captured_output() as (out, err):
                badread.__main__.main(output=self.null)
        out, err = out.getvalue(), err.getvalue()
        self.assertTrue(out.startswith('@'))

    def test_error_model(self):
        test_args = ['badread', 'error_model', '--reference', self.ref_filename,
                     '--reads', self.reads_filename, '--alignment', self.paf_filename]
        with unittest.mock.patch.object(sys, 'argv', test_args):
            with badread.misc.captured_output() as (out, err):
                badread.__main__.main(output=self.null)
        out, err = out.getvalue(), err.getvalue()
        self.assertTrue('AAAAAAT' in out)

    def test_qscore_model(self):
        test_args = ['badread', 'qscore_model', '--reference', self.ref_filename,
                     '--reads', self.reads_filename, '--alignment', self.paf_filename]
        with unittest.mock.patch.object(sys, 'argv', test_args):
            with badread.misc.captured_output() as (out, err):
                badread.__main__.main(output=self.null)
        out, err = out.getvalue(), err.getvalue()
        self.assertTrue('overall;' in out)

    def test_plot(self):
        test_args = ['badread', 'plot', '--reference', self.ref_filename,
                     '--reads', self.reads_filename, '--alignment', self.paf_filename, '--no_plot']
        with unittest.mock.patch.object(sys, 'argv', test_args):
            with badread.misc.captured_output() as (out, err):
                badread.__main__.main()
        out, err = out.getvalue(), err.getvalue()
        self.assertTrue('read_1:' in out)


class TestPythonVersion(unittest.TestCase):

    def test_good_version_1(self):
        with unittest.mock.patch.object(sys, 'version_info') as v_info:
            v_info.major, v_info.minor = 3, 6
            try:
                badread.__main__.check_python_version()
            except SystemExit:
                self.fail('check_python_version() failed when it should not have')

    def test_good_version_2(self):
        with unittest.mock.patch.object(sys, 'version_info') as v_info:
            v_info.major, v_info.minor = 3, 7
            try:
                badread.__main__.check_python_version()
            except SystemExit:
                self.fail('check_python_version() failed when it should not have')

    def test_bad_version_1(self):
        with unittest.mock.patch.object(sys, 'version_info') as v_info:
            v_info.major, v_info.minor = 3, 5
            with self.assertRaises(SystemExit) as cm:
                badread.__main__.check_python_version()
            self.assertNotEqual(cm.exception.code, 0)
            self.assertTrue('Badread requires Python' in str(cm.exception))

    def test_bad_version_2(self):
        with unittest.mock.patch.object(sys, 'version_info') as v_info:
            v_info.major, v_info.minor = 2, 7
            with self.assertRaises(SystemExit) as cm:
                badread.__main__.check_python_version()
            self.assertNotEqual(cm.exception.code, 0)
            self.assertTrue('Badread requires Python' in str(cm.exception))


class TestCommandLine(unittest.TestCase):

    def setUp(self):
        self.ref = os.path.join(os.path.dirname(__file__), 'test_ref_1.fasta')

    def check_simulate_error(self, args, message, ref=None):
        # Ensures that the args throw an error containing the given message.
        if ref is None:
            ref = self.ref
        args = ['simulate', '--reference', ref, '--quantity', '50x'] + args
        args = badread.__main__.parse_args(args)
        with self.assertRaises(SystemExit) as cm:
            badread.__main__.check_simulate_args(args)
        self.assertNotEqual(cm.exception.code, 0)
        self.assertTrue(message in str(cm.exception))

    def check_simulate_no_error(self, args):
        args = ['simulate', '--reference', self.ref, '--quantity', '50x'] + args
        args = badread.__main__.parse_args(args)
        try:
            badread.__main__.check_simulate_args(args)
        except SystemExit:
            self.fail()

    def test_no_args(self):
        # When called with no arguments, the program behaves like it's hit an error: the help text
        # goes to stderr and it returns a non-zero exit code.
        with badread.misc.captured_output() as (out, err):
            with self.assertRaises(SystemExit) as cm:
                badread.__main__.parse_args([])
        self.assertNotEqual(cm.exception.code, 0)
        self.assertTrue('Badread' in err.getvalue().strip())

    def test_help_1(self):
        # When called with -h or --help, the program doesn't behaves like it's hit an error: the
        # help text goes to stdout and it returns a exit code of zero.
        with badread.misc.captured_output() as (out, err):
            with self.assertRaises(SystemExit) as cm:
                badread.__main__.parse_args(['-h'])
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(err.getvalue().strip(), '')
        self.assertTrue('Badread' in out.getvalue().strip())

    def test_help_2(self):
        # Same as above, but with --help.
        with badread.misc.captured_output() as (out, err):
            with self.assertRaises(SystemExit) as cm:
                badread.__main__.parse_args(['--help'])
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(err.getvalue().strip(), '')
        self.assertTrue('Badread' in out.getvalue().strip())

    def test_help_no_colours(self):
        with unittest.mock.patch('subprocess.check_output') as check_output_mock:
            check_output_mock.side_effect = ValueError()
            with badread.misc.captured_output() as (out, err):
                with self.assertRaises(SystemExit) as cm:
                    badread.__main__.parse_args(['--help'])
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(err.getvalue().strip(), '')
        self.assertTrue('Badread' in out.getvalue().strip())

    def test_simulate_help(self):
        with badread.misc.captured_output() as (out, err):
            with self.assertRaises(SystemExit) as cm:
                badread.__main__.parse_args(['simulate', '-h'])
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(err.getvalue().strip(), '')
        self.assertTrue('simulate' in out.getvalue().strip())

    def test_simulate_help_2(self):
        with badread.misc.captured_output() as (out, err):
            with self.assertRaises(SystemExit) as cm:
                badread.__main__.parse_args(['simulate'])
        self.assertNotEqual(cm.exception.code, 0)
        self.assertEqual(out.getvalue().strip(), '')
        self.assertTrue('simulate' in err.getvalue().strip())

    def test_simulate_args_1(self):
        args = badread.__main__.parse_args(['simulate', '--reference', self.ref,
                                           '--quantity', '50x', '--length', '5432,123',
                                           '--identity', '78,89,4', '--glitches', '9876,12,34'])
        self.assertEqual(args.reference, self.ref)
        self.assertEqual(args.quantity, '50x')
        badread.__main__.check_simulate_args(args)
        self.assertEqual(args.mean_frag_length, 5432)
        self.assertEqual(args.frag_length_stdev, 123)
        self.assertEqual(args.mean_identity, 78)
        self.assertEqual(args.max_identity, 89)
        self.assertEqual(args.identity_stdev, 4)
        self.assertEqual(args.glitch_rate, 9876)
        self.assertEqual(args.glitch_size, 12)
        self.assertEqual(args.glitch_skip, 34)

    def test_simulate_args_2(self):
        args = badread.__main__.parse_args(['simulate', '--reference', self.ref,
                                           '--quantity', '250M', '--length', '2345,321',
                                           '--identity', '20,2', '--glitches', '6789,21,43'])
        self.assertEqual(args.reference, self.ref)
        self.assertEqual(args.quantity, '250M')
        badread.__main__.check_simulate_args(args)
        self.assertEqual(args.mean_frag_length, 2345)
        self.assertEqual(args.frag_length_stdev, 321)
        self.assertEqual(args.mean_identity, 20)
        self.assertEqual(args.max_identity, None)
        self.assertEqual(args.identity_stdev, 2)
        self.assertEqual(args.glitch_rate, 6789)
        self.assertEqual(args.glitch_size, 21)
        self.assertEqual(args.glitch_skip, 43)

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
        self.check_simulate_error(['--identity', '20,19,18,17'], 'could not parse')

    def test_bad_simulated_args_12(self):
        self.check_simulate_error(['--identity', '85,101,5'], 'cannot be more than')

    def test_bad_simulated_args_13(self):
        self.check_simulate_error(['--identity', '101,102,5'], 'cannot be more than')

    def test_bad_simulated_args_14(self):
        self.check_simulate_error(['--identity', '85,80,5'], 'cannot be larger than max')

    def test_bad_simulated_args_15(self):
        self.check_simulate_error(['--identity', '85,90,-3'], 'cannot be negative')

    def test_bad_simulated_args_16(self):
        self.check_simulate_error(['--identity', '30,90,5'], 'must be at least')

    def test_bad_simulated_args_17(self):
        self.check_simulate_error(['--identity', '90,30,5'], 'must be at least')

    def test_bad_simulated_args_18(self):
        self.check_simulate_error(['--identity', '1,1'], 'must be at least')

    def test_bad_simulated_args_19(self):
        self.check_simulate_error(['--glitches', 'abc'], 'could not parse')

    def test_bad_simulated_args_20(self):
        self.check_simulate_error(['--glitches', '500,Z,10'], 'could not parse')

    def test_bad_simulated_args_21(self):
        self.check_simulate_error(['--glitches', '500,-10,10'], 'must contain non-negative')

    def test_bad_simulated_args_22(self):
        self.check_simulate_error(['--error_model', 'not_a_file'], 'or a filename')

    def test_bad_simulated_args_23(self):
        self.check_simulate_error(['--qscore_model', 'not_a_file'], 'or a filename')

    def test_bad_simulated_args_24(self):
        self.check_simulate_error([], 'not a file', ref='not_a_file')

    def test_bad_simulated_args_25(self):
        self.check_simulate_error(['--start_adapter_seq', 'ACGTQ'], 'must be a DNA sequence')

    def test_bad_simulated_args_26(self):
        self.check_simulate_error(['--end_adapter_seq', '12ACT'], 'must be a DNA sequence')

    def test_good_simulated_args_1(self):
        self.check_simulate_no_error(['--start_adapter_seq', '12'])

    def test_good_simulated_args_2(self):
        self.check_simulate_no_error(['--end_adapter_seq', '321'])

    def test_good_simulated_args_3(self):
        self.check_simulate_no_error(['--start_adapter_seq', 'ACGCTGCATCG'])

    def test_good_simulated_args_4(self):
        self.check_simulate_no_error(['--end_adapter_seq', 'AAAAAAAAAAAA'])
