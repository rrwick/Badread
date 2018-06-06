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

import edlib
import numpy as np
import random
import sys
import uuid
from .misc import load_fasta, get_random_sequence, reverse_complement, random_chance
from .error_model import ErrorModel, identity_from_edlib_cigar
from .fragment_lengths import FragmentLengths
from .identities import Identities


def simulate(args):
    ref_seqs, ref_depths, ref_circular = load_fasta(args.reference)
    rev_comp_ref_seqs = {name: reverse_complement(seq) for name, seq in ref_seqs.items()}
    target_size = get_target_size(ref_seqs, args.quantity)
    frag_lengths = FragmentLengths(args.lengths, args.mean_frag_length, args.frag_length_stdev)
    error_model = ErrorModel(args.error_model)
    identities = Identities(args.identities, args.mean_identity, args.identity_shape,
                            args.max_identity, error_model)

    start_adapt_rate, start_adapt_amount = adapter_parameters(args.start_adapter_params)
    end_adapt_rate, end_adapt_amount = adapter_parameters(args.end_adapter_params)
    ref_contigs, ref_contig_weights = get_ref_contig_weights(ref_seqs, ref_depths)

    print('', file=sys.stderr)
    print('Generating reads', file=sys.stderr)
    total_size = 0
    while total_size < target_size:
        fragment = [get_start_adapter(start_adapt_rate, start_adapt_amount, args.start_adapter)]
        info = []
        frag_seq, frag_info = get_fragment(frag_lengths, ref_seqs, rev_comp_ref_seqs,
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
            frag_seq, frag_info = get_fragment(frag_lengths, ref_seqs, rev_comp_ref_seqs,
                                               ref_contigs, ref_contig_weights, ref_circular, args)
            fragment.append(frag_seq)
            info += frag_info
        fragment.append(get_end_adapter(end_adapt_rate, end_adapt_amount, args.end_adapter))
        fragment = ''.join(fragment)

        read_identity = identities.get_identity()
        seq, quals = sequence_fragment(fragment, read_identity, args.glitches, args.skips,
                                       error_model)
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


def get_fragment(frag_lengths, ref_seqs, rev_comp_ref_seqs, ref_contigs, ref_contig_weights,
                 ref_circular, args):
    fragment_length = frag_lengths.get_fragment_length()
    fragment_type = get_fragment_type(args)
    if fragment_type == 'junk':
        return get_junk_fragment(fragment_length), ['junk_seq']
    elif fragment_type == 'random':
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


def get_fragment_type(args):
    """
    Returns either 'junk_seq', 'random_seq' or 'good'
    """
    random_draw = random.random()
    if random_draw < args.junk_read_rate:
        return 'junk'
    elif random_draw < args.junk_read_rate + args.random_read_rate:
        return 'random'
    else:
        return 'good'


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


def sequence_fragment(fragment, target_identity, glitches, skips, error_model):

    # TODO: add glitches
    # TODO: add skips

    if error_model.type == 'perfect':
        q_string = ''.join(random.choice('ABCDEFGHI') for _ in range(len(fragment)))
        return fragment, q_string

    # Buffer the fragment a bit so errors can be added to the first and last bases.
    k_size = error_model.kmer_size
    fragment = get_random_sequence(k_size) + fragment + get_random_sequence(k_size)
    frag_len = len(fragment)

    # A list to hold the bases for the errors-added fragment. Note that these values can be ''
    # (meaning the base was deleted) or more than one base (meaning there was an insertion).
    new_fragment_bases = [x for x in fragment]

    errors = 0.0
    change_count = 0

    # # The target identity is the threshold below which we stop adding errors. It is therefore a
    # # little bit higher than the identity we're aiming for.
    # target_identity = 1.0 - ((1.0 - target_identity) * ((frag_len - 1) / frag_len))

    max_kmer_index = len(new_fragment_bases) - 1 - k_size
    while True:

        # To gauge the identity, we first use the number of changes we've added to the fragment,
        # which will probably under-estimate the identity, but it's fast.
        estimated_identity = 1.0 - (errors / frag_len)
        if estimated_identity <= target_identity:
            break

        i = random.randint(0, max_kmer_index)
        kmer = fragment[i:i+k_size]
        new_kmer = error_model.add_errors_to_kmer(kmer)

        # If the error model didn't make any changes (quite common with a non-random error model),
        # we just try again at a different position.
        if kmer == ''.join(new_kmer):
            continue

        for j in range(k_size):
            fragment_base = fragment[i+j]
            new_base = new_kmer[j]  # can actually be more than one base, in cases of insertion

            # If this base is changed in the k-mer and hasn't already been changed, then we apply
            # the change.
            if new_base != fragment_base and fragment_base == new_fragment_bases[i+j]:
                new_fragment_bases[i+j] = new_base
                change_count += 1
                if len(new_base) < 2:  # deletion or substitution
                    errors += 1
                else:  # insertion
                    errors += len(new_base) - 1

                # Every now and then we actually align the new sequence to its original to get a
                # more accurate estimate of edit distance.
                # TODO: this is too slow! I need something more efficient...
                if change_count % 10 == 0:
                    cigar = edlib.align(fragment, ''.join(new_fragment_bases), task='path')['cigar']
                    actual_identity = identity_from_edlib_cigar(cigar)
                    errors = (1.0 - actual_identity) * frag_len

    # Remove the buffer bases.
    new_fragment_bases = new_fragment_bases[k_size:-k_size]
    seq = ''.join(new_fragment_bases)

    qual = 'A' * len(seq)
    # TODO: Phred quality scores
    # TODO: Phred quality scores
    # TODO: Phred quality scores
    # TODO: Phred quality scores
    # TODO: Phred quality scores

    return seq, qual


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
