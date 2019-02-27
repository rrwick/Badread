"""
This module contains functions for drawing gamma and beta distribution histograms to the terminal.

Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Badread

This file is part of Badread. Badread is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Badread is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Badread.
If not, see <http://www.gnu.org/licenses/>.

This is a modified version of the quickhist Python script (github.com/nk412/quickhist).
Original license follows:

Copyright (c) 2015 Nagarjuna Kumarappan

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Author: Nagarjuna Kumarappan <nagarjuna.412@gmail.com>
"""

import math
import numpy as np
import scipy.special
import os
import sys


def draw_hist(y, shape, bins, height, x_tick_interval, y_label='', y_label_space=0,
              print_labels=True, output=sys.stderr):
    max_count = max(y)
    normed_hist_list = [float(x) * height / max_count for x in y]

    # Build plot from top level
    i = 0
    for depth in range(height-1, -1, -1):
        if 0 <= i-y_label_space < len(y_label):
            print(y_label[i-2], end='', file=output)
        else:
            print(' ', end='', file=output)
        print(' │', end='', file=output)

        for item in normed_hist_list:
            floored_item = math.floor(item)
            if floored_item >= depth:
                if floored_item == depth and 0.75 > item % 1 > 0.25:
                    print('\u2596', end='', file=output)  # half bar
                elif floored_item == depth and item % 1 > 0.75:
                    print('\u258c', end='', file=output)  # full bar
                elif floored_item > depth:
                    print('\u258c', end='', file=output)  # full bar
                else:
                    print(' ', end='', file=output)
            else:
                print(' ', end='', file=output)
        print('', file=output)
        i += 1

    # Draw X axis with labels
    line, labels = '  ', '  '
    label = shape[0]
    bin_size = (shape[1] - shape[0]) / bins
    label_step = int(x_tick_interval * bin_size)
    for i in range(bins+1):
        if i == 0:
            line += '├'
            labels += str(label)
        elif i % x_tick_interval == 0:
            line += '┐' if i == bins else '┬'
            label += label_step
            labels += str(label)
        else:
            line += '─'
            labels += ' ' * (len(line) - len(labels))
    print(line, file=output)
    if print_labels:
        print(labels, file=output)


def quickhist_gamma(a, b, n50, height, output=sys.stderr):
    hist_max = int(math.ceil(n50 * 3 / 2000) * 2000)
    tick_interval = 10
    if get_max_width() > 120:
        bin_size = int(hist_max / 100)
    else:
        bin_size = int(hist_max / 50)
    bins = np.asarray([bin_size * (i + 1) for i in range(int(hist_max / bin_size))])
    frags_y, bases_y = [], []
    for b_val in bins:
        x = b_val - (bin_size / 2)  # take the value from the middle of the bin

        # These are the original functions to get the density:
        # frag_y = (b ** a) * (x ** (a-1)) * (np.exp(-x * b) / (scipy.special.gamma(a)))
        # base_y = (b ** (a+1)) * (x ** a) * (np.exp(-x * b) / (scipy.special.gamma(a+1)))
        # But to avoid overflows, I had to log-ify them:
        frag_y = np.exp((-x*b) + ((a-1)*np.log(x)) + (a*np.log(b)) - scipy.special.gammaln(a))
        base_y = np.exp((-x*b) + (a*np.log(x)) + ((a+1)*np.log(b)) - scipy.special.gammaln(a+1))

        frags_y.append(frag_y)
        bases_y.append(base_y)

    shape = (0, hist_max)
    draw_hist(frags_y, shape, len(bins), height, tick_interval, 'frags', 2, print_labels=False,
              output=output)
    draw_hist(bases_y, shape, len(bins), height, tick_interval, 'bases', 2, output=output)


def quickhist_beta(a, b, max_identity, height, output=sys.stderr):
    hist_min, hist_max = 50, 100
    tick_interval = 10
    if get_max_width() > 120:
        bin_size = 0.5
    else:
        bin_size = 1
    bins = (np.arange(hist_min, hist_max, bin_size) + bin_size) / 100 / max_identity
    y = []
    for b_val in bins:
        x = b_val - (bin_size / 200)  # take the value from the middle of the bin
        if x < 1:
            # This is the original function to get the density:
            # beta = x**(a-1) * (1-x)**(b-1) / scipy.special.beta(a, b)
            # But to avoid overflows, I had to log-ify it:
            beta = np.exp((a-1)*np.log(x) + (b-1)*np.log(1-x) - scipy.special.betaln(a, b))
        else:
            beta = 0.0
        y.append(beta)
    bins *= max_identity
    shape = (hist_min, hist_max)
    draw_hist(y, shape, len(bins), height, tick_interval, output=output)


def get_terminal_size_stderr(fallback=(80, 24)):
    """
    Unlike shutil.get_terminal_size, which looks at stdout, this looks at stderr.
    """
    try:
        size = os.get_terminal_size(sys.__stderr__.fileno())
    except (AttributeError, ValueError, OSError):
        size = os.terminal_size(fallback)
    return size


def get_max_width():
    col_count = get_terminal_size_stderr().columns
    if col_count < 80:
        return 80
    if col_count > 160:
        return 160
    return col_count
