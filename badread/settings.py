"""
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Badread

The version is stored here in a separate file so it can exist in only one place.
http://stackoverflow.com/questions/458550

This file is part of Badread. Badread is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Badread is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Badread.
If not, see <http://www.gnu.org/licenses/>.
"""


# When adding errors to a sequence (sequence_fragment function in simulate.py), it is necessary to
# keep a running estimate of the sequence identity. This estimate can be greatly improved by
# actually aligning the error-added sequence to the original sequence, but this is computationally
# expensive. These settings control how often such alignment takes place and how much sequence is
# used in the alignments. Set ALIGNMENT_INTERVAL to a very large value to turn alignments off, so
# the identity estimates will be based on error count alone.
ALIGNMENT_INTERVAL = 25
ALIGNMENT_SIZE = 1000


# I don't let users set a very small minimum mean read length (e.g. 2) or very low minimum read
# identity (e.g. 50%) as that might break some things. These settings control how low they can go.
MIN_MEAN_READ_LENGTH = 100
MIN_MEAN_READ_IDENTITY = 50


# If a random qscore model is used, this controls the range of the qualities.
RANDOM_QSCORE_MAX = 20
