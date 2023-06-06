"""
This module contains a class for describing read identity distributions and related functions.

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

import numpy as np
import sys
from .quickhist import quickhist_beta
from .misc import float_to_str, print_in_two_columns


class Identities(object):

    def __init__(self, mean, stdev, max_identity, output=sys.stderr):
        self.mean, self.stdev, self.max_identity = None, None, None
        self.beta_a, self.beta_b = None, None

        print('', file=output)

        # There are two possible types of identity distributions: a three-parameter beta
        # distribution that describes read identities, or a two-parameter normal distribution that
        # describes read qscores.
        if max_identity is None:
            self.type = "normal"
            self.set_up_normal(mean, stdev, output)
        else:
            self.type = "beta"
            self.set_up_beta(mean, stdev, max_identity, output)

    def set_up_beta(self, mean, stdev, max_identity, output):
        # Divide by 100 to convert from percentage to fraction
        self.mean = mean / 100.0
        self.stdev = stdev / 100.0
        self.max_identity = max_identity / 100.0

        if self.mean == self.max_identity:
            self.beta_a, self.beta_b = None, None
            print(f'Using a constant read identity of {self.mean * 100}%', file=output)
        elif self.stdev == 0.0:
            self.max_identity = self.mean
            print(f'Using a constant read identity of {self.mean * 100}%', file=output)
        else:
            print('Generating read identities from a beta distribution:', file=output)
            self.beta_a, self.beta_b = beta_parameters(mean, stdev, max_identity)
            print_in_two_columns(f'  mean  = {float_to_str(self.mean * 100):>3}%',
                                 f'  max   = {float_to_str(self.max_identity * 100):>3}%',
                                 f'  stdev = {float_to_str(self.stdev * 100):>3}%',
                                 'shape parameters:',
                                 f'  alpha = {self.beta_a:.4e}',
                                 f'  beta  = {self.beta_b:.4e}',
                                 output=output)
            quickhist_beta(self.beta_a, self.beta_b, self.max_identity, 8, output=output)

    def set_up_normal(self, mean, stdev, output):
        self.mean = mean
        self.stdev = stdev

        if self.stdev == 0.0:
            self.max_identity = self.mean
            print(f'Using a constant read qscore of {self.mean}', file=output)
        else:
            print('Generating read qscores from a normal distribution:', file=output)
            print(f'  mean  = {float_to_str(self.mean):>3}', file=output)
            print(f'  stdev = {float_to_str(self.stdev):>3}', file=output)

    def get_identity(self):
        while True:
            if self.type == "beta":
                identity = self.get_beta_identity()
            else:
                identity = self.get_normal_identity()
            if 0 <= identity <= 100:
                return identity

    def get_beta_identity(self):
        if self.mean == self.max_identity:
            return self.mean
        else:  # beta distribution
            return self.max_identity * np.random.beta(self.beta_a, self.beta_b)

    def get_normal_identity(self):
        qscore = np.random.normal(self.mean, self.stdev)
        return 1.0 - 10**(-qscore / 10)


def beta_parameters(beta_mean, beta_stdev, beta_max):
    u, s, m = beta_mean, beta_stdev, beta_max
    beta_a = (((1-(u/m)) / ((s/m)**2)) - (m/u)) * ((u/m)**2)
    beta_b = beta_a * ((m/u) - 1)
    if beta_a < 0.0 or beta_b < 0.0:
        sys.exit('Error: invalid beta parameters for identity distribution - trying increasing '
                 'the maximum identity or reducing the standard deviation')
    return beta_a, beta_b
