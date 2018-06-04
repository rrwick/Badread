"""
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Badread

This file is part of Badread. Badread is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Badread is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Badread.
If not, see <http://www.gnu.org/licenses/>.

This is a modified version of the quickhist Python script (github.com/nk412/quickhist). Original
license follows:

Copyright (c) 2015 Nagarjuna Kumarappan

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Author: Nagarjuna Kumarappan <nagarjuna.412@gmail.com>
"""


import math
import numpy as np
import scipy.special
import sys

# Constants
CHARS = {'full_top': '\u258c', 'half_top': '\u2596', 'fill': '\u258c'}


def draw_hist(normed_hist_list, shape, bins, height):

    # TODO: add x axis ticks so you can see exactly where the numbers are
    # TODO: add more x axis number labels

    """ takes a list of histogram bin counts and shape(min,max) of input """

    # Build plot from top level
    for depth in range(height-1, -1, -1):

        # Draw Y axis
        sys.stderr.write(' │')

        # Draw bars
        for item in normed_hist_list:
            floored_item = math.floor(item)
            if floored_item >= depth:
                if floored_item == depth and 0.75 > item % 1 > 0.25:
                    sys.stderr.write(CHARS['half_top'])
                elif floored_item == depth and item % 1 > 0.75:
                    sys.stderr.write(CHARS['full_top'])
                elif floored_item > depth:
                    sys.stderr.write(CHARS['fill'])
                else:
                    sys.stderr.write(' ')
                continue
            else:
                sys.stderr.write(' ')
        print('', file=sys.stderr)

    # Draw X axis 
    print(' └' + '─'*(bins+2), file=sys.stderr)
    print(' ' + str(shape[0]) + ' '*(bins-3) + str(shape[1]), file=sys.stderr)


def quickhist(bins, height, hist_min, hist_max, input_list):
    hist_list, bin_edges = np.histogram(input_list, range=(hist_min, hist_max), bins=bins)
    max_count = max(hist_list)
    shape = (hist_min, hist_max)
    normed_hist_list = [float(x)*height/max_count for x in hist_list]
    
    # Draw the histogram
    draw_hist(normed_hist_list, shape, bins, height)


def quickhist_gamma(k, t, n50, height):
    width = get_max_width()
    bin_size = 250
    while True:
        hist_max = int(math.ceil(n50 * 2 / bin_size) * bin_size)
        if hist_max / bin_size <= width:
            break
        bin_size *= 2

    bins = np.asarray([bin_size * (i + 1) for i in range(int(hist_max / bin_size))])
    y = bins ** (k - 1) * (np.exp(-bins / t) / (scipy.special.gamma(k) * t ** k))

    # TODO: tweak things such that the values plotted are from the middle of the bins.

    max_count = max(y)
    shape = (0, hist_max)
    normed_hist_list = [float(x) * height / max_count for x in y]
    draw_hist(normed_hist_list, shape, len(bins), height)


def quickhist_beta(a, b, max_identity, height):
    width = get_max_width()
    hist_min, hist_max = 50, 100
    bin_size = 0.5
    while True:
        if (hist_max - hist_min) / bin_size <= width:
            break
        bin_size *= 2

    bins = (np.arange(hist_min, hist_max, bin_size) + bin_size) / 100 / max_identity
    y = []
    for x in bins:
        if x < 1:
            y.append(x**(a-1) * (1-x)**(b-1) / scipy.special.beta(a, b))
        else:
            y.append(0.0)
    bins *= max_identity

    # TODO: tweak things such that the values plotted are from the middle of the bins.

    max_count = max(y)
    shape = (hist_min, hist_max)
    normed_hist_list = [float(x) * height / max_count for x in y]
    draw_hist(normed_hist_list, shape, len(bins), height)


def get_max_width():
    # TODO: automatically set width, based on the terminal width
    return 70
