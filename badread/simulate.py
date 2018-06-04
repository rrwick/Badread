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
import random
import scipy.special
import sys
import uuid
from .misc import load_fasta, get_random_sequence, reverse_complement, random_chance
from .error_model import ErrorModel


def simulate(args):
    ref_seqs, ref_depths, ref_circular = load_fasta(args.reference)
    rev_comp_ref_seqs = {name: reverse_complement(seq) for name, seq in ref_seqs.items()}
    target_size = get_target_size(ref_seqs, args.quantity)
    gamma_k, gamma_t = gamma_parameters(args.mean_read_length, args.read_length_stdev)
    beta_a, beta_b = beta_parameters(args.mean_read_identity, args.read_identity_dist_shape,
                                     args.max_read_identity)

    start_adapt_rate, start_adapt_amount = adapter_parameters(args.start_adapter_params)
    end_adapt_rate, end_adapt_amount = adapter_parameters(args.end_adapter_params)
    ref_contigs, ref_contig_weights = get_ref_contig_weights(ref_seqs, ref_depths)

    error_model = ErrorModel(args.error_model)

    total_size = 0
    while total_size < target_size:
        fragment = [get_start_adapter(start_adapt_rate, start_adapt_amount, args.start_adapter)]
        info = []
        frag_seq, frag_info = get_fragment(gamma_k, gamma_t, ref_seqs, rev_comp_ref_seqs,
                                           ref_contigs, ref_contig_weights, ref_circular, args)
        fragment.append(frag_seq)
        info += frag_info

        while random_chance(args.chimera_rate):
            info.append('chimera')
            if random_chance(0.25):
                fragment.append(get_end_adapter(end_adapt_rate, end_adapt_amount, args.end_adapter))
            if random_chance(0.25):
                fragment.append(get_start_adapter(start_adapt_rate, start_adapt_amount,
                                                  args.start_adapter))
            frag_seq, frag_info = get_fragment(gamma_k, gamma_t, ref_seqs, rev_comp_ref_seqs,
                                               ref_contigs, ref_contig_weights, ref_circular, args)
            fragment.append(frag_seq)
            info += frag_info
        fragment.append(get_end_adapter(end_adapt_rate, end_adapt_amount, args.end_adapter))
        fragment = ''.join(fragment)

        seq, quals = sequence_fragment(fragment, beta_a, beta_b, args.max_read_identity,
                                       args.glitches, args.skips, error_model)
        read_name = uuid.uuid4()

        print('@{} {}'.format(read_name, ','.join(info)))
        print(seq)
        print('+')
        print(quals)

        total_size += len(fragment)

    # print('', file=sys.stderr)
    # print('N50: {}'.format(get_n50(reads)), file=sys.stderr)
    # print('Mean identity: {:.4f}'.format(get_mean_id(reads)), file=sys.stderr)


def get_ref_contig_weights(ref_seqs, ref_depths):
    ref_contigs = [x[0] for x in ref_depths]
    ref_contig_weights = [x[1] * len(ref_seqs[x[0]]) for x in ref_depths]
    return ref_contigs, ref_contig_weights


def get_target_size(reference, quantity):
    try:
        return int(quantity)
    except ValueError:
        pass
    quantity = quantity.lower()
    try:
        last_char = quantity[-1]
        value = float(quantity[:-1])
        if last_char == 'x':
            return int(round(value * sum(len(x) for x in reference.values())))
        elif last_char == 'g':
            return int(round(value * 1000000000))
        elif last_char == 'm':
            return int(round(value * 1000000))
        elif last_char == 'k':
            return int(round(value * 1000))
    except (ValueError, IndexError):
        pass
    sys.exit('Error: could not parse quantity\n'
             '--quantity must be either an absolute value (e.g. 250M) or a relative depth '
             '(e.g. 25x)')


def get_fragment(gamma_k, gamma_t, ref_seqs, rev_comp_ref_seqs, ref_contigs, ref_contig_weights,
                 ref_circular, args):
    fragment_length = int(round(np.random.gamma(gamma_k, gamma_t)))
    random_draw = random.random()
    if random_draw < args.junk_read_rate:
        return get_junk_fragment(fragment_length), ['junk_seq']
    elif random_draw < args.junk_read_rate + args.random_read_rate:
        return get_random_sequence(fragment_length), ['random_seq']

    # The get_real_fragment function can return nothing (due to --lose_small_plasmids) so we try
    # repeatedly until we get a result.
    for _ in range(100):
        seq, info = get_real_fragment(fragment_length, ref_seqs, rev_comp_ref_seqs, ref_contigs,
                                      ref_contig_weights, ref_circular, args)
        if seq != '':
            return seq, info
    sys.exit('Error: failed to generate any sequence fragments - are your read lengths '
             'incompatible with your reference contig lengths?')


def get_real_fragment(fragment_length, ref_seqs, rev_comp_ref_seqs, ref_contigs,
                      ref_contig_weights, ref_circular, args):

    if len(ref_contigs) == 1:
        contig = ref_contigs[0]
    else:
        contig = random.choices(ref_contigs, weights=ref_contig_weights)[0]
    info = [contig]
    if random_chance(0.5):
        seq = ref_seqs[contig]
        info.append('+strand')
    else:
        seq = rev_comp_ref_seqs[contig]
        info.append('-strand')

    # If the reference contig is linear and the fragment length is long enough, then we just
    # return the entire fragment, start to end.
    if fragment_length >= len(seq) and not ref_circular[contig]:
        info.append('0-' + str(len(seq)))
        return seq

    # If the reference contig is circular and the fragment length is too long, then we either
    # fail to get the read (if --lose_small_plasmids was used) or bring the fragment size back
    # down to the contig size.
    if fragment_length > len(seq) and ref_circular[contig]:
        if args.lose_small_plasmids:
            return '', ''
        else:
            fragment_length = len(seq)

    start_pos = random.randint(0, len(seq)-1)
    end_pos = start_pos + fragment_length

    info.append('{}-{}'.format(start_pos, end_pos))

    # For circular contigs, we may have to loop the read around the contig.
    if ref_circular[contig]:
        if end_pos <= len(seq):
            return seq[start_pos:end_pos], info
        else:
            looped_end_pos = end_pos - len(seq)
            assert looped_end_pos > 0
        return seq[start_pos:] + seq[:looped_end_pos], info

    # For linear contigs, we don't care if the ending position is off the end - that will just
    # result in the read ending at the sequence end (and being shorter than the fragment
    # length).
    else:
        return seq[start_pos:end_pos], info


def get_junk_fragment(fragment_length):
    repeat_length = random.randint(1, 8)
    repeat_count = int(round(fragment_length / repeat_length))
    return get_random_sequence(repeat_length) * repeat_count


def sequence_fragment(fragment, beta_a, beta_b, beta_max, glitches, skips, error_model):

    # TODO: add glitches
    # TODO: add skips

    if error_model.type == 'perfect':
        q_string = ''.join(random.choice('ABCDEFGHI') for _ in range(len(fragment)))
        return fragment, q_string

    read_identity = beta_max * np.random.beta(beta_a, beta_b)
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO
    return fragment, 'A' * len(fragment)


def get_start_adapter(rate, amount, adapter):
    if random_chance(rate):
        if amount == 1.0:
            return adapter
        adapter_frag_length = get_adapter_frag_length(amount, adapter)
        start_pos = len(adapter) - adapter_frag_length
        return adapter[start_pos:]
    return ''


def get_end_adapter(rate, amount, adapter):
    if random_chance(rate):
        if amount == 1.0:
            return adapter
        adapter_frag_length = get_adapter_frag_length(amount, adapter)
        return adapter[:adapter_frag_length]
    return ''


def get_adapter_frag_length(amount, adapter):
    """
    Uses a beta distribution to choose the fragment length.
    desmos.com/calculator/l2ssxwgqqf
    """
    beta_a = 2.0 * amount
    beta_b = 2.0 - beta_a
    return round(int(len(adapter) * np.random.beta(beta_a, beta_b)))


def gamma_parameters(gamma_mean, gamma_stdev):
    # Shape and rate parametrisation:
    gamma_a = (gamma_mean ** 2) / (gamma_stdev ** 2)
    gamma_b = gamma_mean / (gamma_stdev ** 2)

    # Shape and scale parametrisation:
    gamma_k = (gamma_mean ** 2) / (gamma_stdev ** 2)
    gamma_t = (gamma_stdev ** 2) / gamma_mean

    print('', file=sys.stderr)
    print('Read length gamma distribution parametrisation:', file=sys.stderr)
    print('  shape: alpha = k = ' + '%.4e' % gamma_a, file=sys.stderr)
    print('  rate:   beta = ' + '%.4e' % gamma_b, file=sys.stderr)
    print('  scale: theta = ' + '%.4e' % gamma_t, file=sys.stderr)

    print('', file=sys.stderr)
    print('Theoretical read lengths:', file=sys.stderr)
    print('  N10: {} bp'.format(int(round(find_n_value(gamma_a, gamma_b, 10)))), file=sys.stderr)
    print('  N50: {} bp'.format(int(round(find_n_value(gamma_a, gamma_b, 50)))), file=sys.stderr)
    print('  N90: {} bp'.format(int(round(find_n_value(gamma_a, gamma_b, 90)))), file=sys.stderr)

    return gamma_k, gamma_t


def beta_parameters(beta_mean, beta_shape, beta_max):
    beta_a = (beta_mean / beta_max) * (beta_shape**2)
    beta_b = (1 - (beta_mean / beta_max)) * (beta_shape**2)

    print('', file=sys.stderr)
    print('Read identity beta distribution parametrisation:', file=sys.stderr)
    print('  a = ' + '%.4e' % beta_a, file=sys.stderr)
    print('  b = ' + '%.4e' % beta_b, file=sys.stderr)

    return beta_a, beta_b


def find_n_value(a, b, n):
    target = 1.0 - (n / 100.0)
    bottom_range = 0.0
    top_range = 1.0
    while base_distribution_integral(a, b, top_range) < target:
        bottom_range = top_range
        top_range *= 2
    guess = (bottom_range + top_range) / 2.0
    while True:
        integral = base_distribution_integral(a, b, guess)
        if top_range - bottom_range < 0.01:
            return guess
        if guess == target:
            return guess
        elif integral < target:
            bottom_range = guess
            guess = (bottom_range + top_range) / 2.0
        else:  # integral > target:
            top_range = guess
            guess = (bottom_range + top_range) / 2.0


def base_distribution_integral(a, b, x):
    # TODO: this function bombs out if the value of a is too large.
    # Could I use log-gamma functions to avoid this, perhaps?
    g = scipy.special.gamma(a+1)
    h = inc_gamma(a+1, b*x)
    return (g - h) / g


def inc_gamma(a, b):
    """
    SciPy seems to define the incomplete Gamma function a bit differently than WolframAlpha (which
    I used to do the calc), so this function should represent a WolframAlpha incomplete Gamma.
    https://stackoverflow.com/questions/38713199/incomplete-gamma-function-in-scipy
    """
    return scipy.special.gamma(a) * (1-scipy.special.gammainc(a, b))


def get_n50(reads):
    read_lengths = sorted([x[0] for x in reads], reverse=True)
    n50_target = sum(read_lengths) / 2
    n50, sum_so_far = 0, 0
    for read_length in read_lengths:
        sum_so_far += read_length
        if sum_so_far >= n50_target:
            return read_length


def get_mean_id(reads):
    total_length, total_id = 0, 0.0
    for read_length, read_id in reads:
        total_length += read_length
        total_id += (read_length * read_id)
    return total_id / total_length


def adapter_parameters(param_str):
    parts = param_str.split(',')
    if len(parts) == 2:
        try:
            return [float(x) for x in parts]
        except ValueError:
            pass
    sys.exit('Error: adapter parameters must be two comma-separated values between 0 and 1')
