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

import statistics
import unittest

import badread.simulate
import badread.identities
import badread.error_model
import badread.misc


class TestGlitches(unittest.TestCase):

    def setUp(self):
        self.frag_length = 1000
        self.trials = 100

    def test_no_glitches_1(self):
        for i in range(self.trials):
            frag = badread.misc.get_random_sequence(self.frag_length)
            glitched_frag = badread.simulate.add_glitches(frag, 0, 0, 0)
            self.assertEqual(frag, glitched_frag)

    def test_no_glitches_2(self):
        # Having a positive glitch rate shouldn't matter if the size and skip are 0.
        for i in range(self.trials):
            frag = badread.misc.get_random_sequence(self.frag_length)
            glitched_frag = badread.simulate.add_glitches(frag, 100, 0, 0)
            self.assertEqual(frag, glitched_frag)

    def test_no_skip_glitches(self):
        # If the glitches have no skip, then the resulting sequence can only get longer.
        new_lengths = []
        for i in range(self.trials):
            frag = badread.misc.get_random_sequence(self.frag_length)
            glitched_frag = badread.simulate.add_glitches(frag, 100, 10, 0)
            self.assertTrue(len(glitched_frag) >= self.frag_length)
            new_lengths.append(len(glitched_frag))
        self.assertTrue(statistics.mean(new_lengths) > self.frag_length)

    def test_no_size_glitches(self):
        # If the glitches have no size, then the resulting sequence can only get shorter.
        new_lengths = []
        for i in range(self.trials):
            frag = badread.misc.get_random_sequence(self.frag_length)
            glitched_frag = badread.simulate.add_glitches(frag, 100, 0, 10)
            self.assertTrue(len(glitched_frag) <= self.frag_length)
            new_lengths.append(len(glitched_frag))
        self.assertTrue(statistics.mean(new_lengths) < self.frag_length)

    def test_size_change(self):
        # When the glitches have size and skip, we expect some fragments to get longer and others
        # to get shorter.
        longer_count, shorter_count = 0, 0
        for i in range(self.trials):
            frag = badread.misc.get_random_sequence(1000)
            glitched_frag = badread.simulate.add_glitches(frag, 100, 10, 10)
            if len(glitched_frag) > self.frag_length:
                longer_count += 1
            elif len(glitched_frag) < self.frag_length:
                shorter_count += 1
        self.assertGreater(longer_count, 0)
        self.assertGreater(shorter_count, 0)

    def test_occurrences(self):
        # When the glitches have a rate the same as the read length, we expect some reads to have
        # glitches and other to not.
        glitch_count, no_glitch_count = 0, 0
        for i in range(self.trials):
            frag = badread.misc.get_random_sequence(1000)
            glitched_frag = badread.simulate.add_glitches(frag, 1000, 10, 10)
            if frag == glitched_frag:
                no_glitch_count += 1
            else:
                glitch_count += 1
        self.assertGreater(glitch_count, 0)
        self.assertGreater(no_glitch_count, 0)

    def test_less_than_one_1(self):
        # Giving a glitch rate/size/skip that's between 0 and 1 used to cause a crash until I fixed
        # the bug.
        for i in range(self.trials):
            frag = badread.misc.get_random_sequence(self.frag_length)
            _ = badread.simulate.add_glitches(frag, 0.5, 10, 10)

    def test_less_than_one_2(self):
        for i in range(self.trials):
            frag = badread.misc.get_random_sequence(self.frag_length)
            _ = badread.simulate.add_glitches(frag, 1000, 0.5, 10)

    def test_less_than_one_3(self):
        for i in range(self.trials):
            frag = badread.misc.get_random_sequence(self.frag_length)
            _ = badread.simulate.add_glitches(frag, 1000, 10, 0.5)
