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
import badread.simulate


class TestPerfectSequenceFragment(unittest.TestCase):
    """
    Tests the sequence_fragment function with a perfect error model.
    """
    def setUp(self):
        self.model = badread.error_model.ErrorModel('perfect')

    def test_sequence_fragment_1(self):
        frag = 'GACCCAGTTTTTTTACTGATTCAGCGTAGGTGCTCTGATCTTCACGCATCTTTGACCGCC'
        seq, qual = badread.simulate.sequence_fragment(frag, 0.0, 0.0, 0.0, 0.0, 0.0, self.model)
        self.assertEqual(frag, seq)
        self.assertEqual(len(frag), len(qual))
        for q in qual:
            self.assertTrue(q in 'ABCDEFGHI')

    def test_sequence_fragment_2(self):
        # Beta distribution parameters are ignored for perfect error models.
        frag = 'TATAAAGACCCCACTTTTGAAGCCAGAGGTAATGGCCGTGATGGCGTTAAATTCCCTTCC'
        seq, qual = badread.simulate.sequence_fragment(frag, 0.5, 0.5, 0.5, 0.0, 0.0, self.model)
        self.assertEqual(frag, seq)
        self.assertEqual(len(frag), len(qual))
        for q in qual:
            self.assertTrue(q in 'ABCDEFGHI')
