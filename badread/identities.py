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
"""

import numpy as np
import sys
from .quickhist import quickhist_beta
from .misc import float_to_str


class Identities(object):

    def __init__(self, mean, shape, max_identity, error_model=None,
                 output=sys.stderr):
        self.mean = mean / 100.0  # convert from percentage
        self.shape = shape
        self.max_identity = max_identity / 100.0  # convert from percentage
        print('', file=output)

        if self.mean == self.max_identity:
            self.beta_a, self.beta_b = None, None
            print('Using a constant read identity of {}%'.format(self.mean * 100), file=output)
        else:  # beta distribution
            print('Generating read identities from a beta distribution:', file=output)
            self.beta_a, self.beta_b = beta_parameters(mean, shape, max_identity)
            print('  alpha (shape) = ' + '%.4e' % self.beta_a, file=output)
            print('  beta (shape)  = ' + '%.4e' % self.beta_b, file=output)
            print('  mean: {}%'.format(float_to_str(self.mean * 100)), file=output)
            print('  max:  {}%'.format(float_to_str(self.max_identity * 100)), file=output)
            print('  shape: {}'.format(float_to_str(self.shape)), file=output)
            quickhist_beta(self.beta_a, self.beta_b, self.max_identity, 8, output=output)

    def get_identity(self):
        if self.mean == self.max_identity:
            return self.mean
        else:  # beta distribution
            return self.max_identity * np.random.beta(self.beta_a, self.beta_b)


def beta_parameters(beta_mean, beta_shape, beta_max):
    beta_a = (beta_mean / beta_max) * (beta_shape**2)
    beta_b = (1 - (beta_mean / beta_max)) * (beta_shape**2)
    return beta_a, beta_b
