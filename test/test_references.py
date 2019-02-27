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

import badread.simulate


class TestLinearFragments(unittest.TestCase):

    def setUp(self):
        null = open(os.devnull, 'w')
        ref_filename = os.path.join(os.path.dirname(__file__), 'test_ref_1.fasta')
        self.ref_seqs, self.ref_depths, self.ref_circular = \
            badread.simulate.load_reference(ref_filename, output=null)
        null.close()

    def test_names(self):
        self.assertEqual(sorted(self.ref_seqs.keys()),
                         ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K'])

    def test_seqs(self):
        self.assertEqual(self.ref_seqs['A'], 'GTCCGCTCAACCTCCTGTTTATCTTAAACC')
        self.assertEqual(self.ref_seqs['E'], 'ATAAGTCCCGCCCTTCGCCCGTTTTCAGGG')
        self.assertEqual(self.ref_seqs['J'], 'AACACACCTTTTTCGAGTCACCCGACGCGC')

    def test_circularity(self):
        self.assertEqual(self.ref_circular['A'], True)
        self.assertEqual(self.ref_circular['B'], False)
        self.assertEqual(self.ref_circular['C'], False)
        self.assertEqual(self.ref_circular['D'], True)
        self.assertEqual(self.ref_circular['E'], False)
        self.assertEqual(self.ref_circular['F'], True)
        self.assertEqual(self.ref_circular['G'], False)
        self.assertEqual(self.ref_circular['H'], False)
        self.assertEqual(self.ref_circular['I'], False)
        self.assertEqual(self.ref_circular['J'], False)
        self.assertEqual(self.ref_circular['K'], False)

    def test_depth(self):
        self.assertEqual(self.ref_depths['A'], 1.0)
        self.assertEqual(self.ref_depths['B'], 1.0)
        self.assertEqual(self.ref_depths['C'], 1.0)
        self.assertEqual(self.ref_depths['D'], 2.0)
        self.assertEqual(self.ref_depths['E'], 3.0)
        self.assertEqual(self.ref_depths['F'], 4.0)
        self.assertEqual(self.ref_depths['G'], 5.0)
        self.assertEqual(self.ref_depths['H'], 6.0)
        self.assertEqual(self.ref_depths['I'], 1.0)
        self.assertEqual(self.ref_depths['J'], 5.4321)
        self.assertEqual(self.ref_depths['K'], 1.23456)
