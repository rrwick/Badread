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


class Identities(object):

    def __init__(self, distribution, mean, shape=None, max_identity=100, error_model=None,
                 output=sys.stderr):
        self.distribution = distribution
        self.mean = mean / 100.0  # convert from percentage
        self.shape = shape
        self.max_identity = max_identity / 100.0  # convert from percentage
        assert distribution == 'constant' or distribution == 'beta'
        print('', file=output)

        if error_model is not None and error_model.type == 'perfect':
            if self.max_identity < 1:
                print('Setting max identity to 100% based on perfect error model', file=output)
                self.max_identity = 1
            if self.mean < 1:
                print('Setting mean identity to 100% based on perfect error model', file=output)
                self.mean = 1
            if self.distribution != 'constant':
                print('Switching to "constant" identity based on perfect error model', file=output)
                self.distribution = 'constant'

        if self.distribution == 'beta' and self.mean == self.max_identity:
            self.distribution = 'constant'
            print('Switching from "beta" to "constant" identity distribution because '
                  'mean identity equals max identity', file=output)
        if self.distribution == 'constant':
            self.beta_a, self.beta_b = None, None
            print('Using a constant read identity of {}%'.format(self.mean * 100), file=output)
        else:  # beta distribution
            print('Generating read identities from a beta distribution:', file=output)
            self.beta_a, self.beta_b = beta_parameters(mean, shape, max_identity)
            print('  alpha (shape) = ' + '%.4e' % self.beta_a, file=output)
            print('  beta (shape)  = ' + '%.4e' % self.beta_b, file=output)
            print('  mean: {}%'.format(self.mean * 100), file=output)
            print('  max:  {}%'.format(self.max_identity * 100), file=output)
            print('  shape: {}'.format(self.shape), file=output)

    def get_identity(self):
        if self.distribution == 'constant':
            return self.mean
        else:  # beta distribution
            return self.max_identity * np.random.beta(self.beta_a, self.beta_b)


def beta_parameters(beta_mean, beta_shape, beta_max):
    beta_a = (beta_mean / beta_max) * (beta_shape**2)
    beta_b = (1 - (beta_mean / beta_max)) * (beta_shape**2)
    return beta_a, beta_b
