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
from .misc import load_fasta, get_random_sequence, reverse_complement, random_chance, \
    float_to_str
from .error_model import ErrorModel, identity_from_edlib_cigar
from .qscore_model import QScoreModel, align_sequences_from_edlib_cigar
from .fragment_lengths import FragmentLengths
from .identities import Identities
from . import settings


def simulate(args):
    ref_seqs, ref_depths, ref_circular = load_fasta(args.reference)
    rev_comp_ref_seqs = {name: reverse_complement(seq) for name, seq in ref_seqs.items()}
    target_size = get_target_size(ref_seqs, args.quantity)
    frag_lengths = FragmentLengths(args.mean_frag_length, args.frag_length_stdev)
    error_model = ErrorModel(args.error_model)
    qscore_model = QScoreModel(args.qscore_model)
    identities = Identities(args.mean_identity, args.identity_shape, args.max_identity, error_model)

    start_adapt_rate, start_adapt_amount = adapter_parameters(args.start_adapter)
    end_adapt_rate, end_adapt_amount = adapter_parameters(args.end_adapter)
    ref_contigs, ref_contig_weights = get_ref_contig_weights(ref_seqs, ref_depths)
    chimera_rate = args.chimeras / 100  # percentage to fraction
    print_glitch_summary(args.glitch_rate, args.glitch_size, args.glitch_skip)

    print('', file=sys.stderr)
    print_progress(0)
    total_size = 0
    while total_size < target_size:
        fragment = [get_start_adapter(start_adapt_rate, start_adapt_amount, args.start_adapter_seq)]
        info = []
        frag_seq, frag_info = get_fragment(frag_lengths, ref_seqs, rev_comp_ref_seqs,
                                           ref_contigs, ref_contig_weights, ref_circular, args)
        fragment.append(frag_seq)
        info += frag_info

        while random_chance(chimera_rate):
            info.append('chimera')
            if random_chance(0.25):
                fragment.append(get_end_adapter(end_adapt_rate, end_adapt_amount,
                                                args.end_adapter_seq))
            if random_chance(0.25):
                fragment.append(get_start_adapter(start_adapt_rate, start_adapt_amount,
                                                  args.start_adapter_seq))
            frag_seq, frag_info = get_fragment(frag_lengths, ref_seqs, rev_comp_ref_seqs,
                                               ref_contigs, ref_contig_weights, ref_circular, args)
            fragment.append(frag_seq)
            info += frag_info
        fragment.append(get_end_adapter(end_adapt_rate, end_adapt_amount, args.end_adapter_seq))
        fragment = ''.join(fragment)
        fragment = add_glitches(fragment, args.glitch_rate, args.glitch_size, args.glitch_skip)
        read_identity = identities.get_identity()
        seq, quals = sequence_fragment(fragment, read_identity, error_model, qscore_model)
        read_name = uuid.uuid4()

        print('@{} {}'.format(read_name, ','.join(info)))
        print(seq)
        print('+')
        print(quals)

        total_size += len(fragment)
        print_progress(total_size)

    print('\n', file=sys.stderr)


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

    # The get_real_fragment function can return nothing (due to --small_plasmid_bias) so we try
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
    junk_read_rate = args.junk_reads / 100      # percentage to fraction
    random_read_rate = args.random_reads / 100  # percentage to fraction
    random_draw = random.random()
    if random_draw < junk_read_rate:
        return 'junk'
    elif random_draw < junk_read_rate + random_read_rate:
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
    # fail to get the read (if --small_plasmid_bias was used) or bring the fragment size back
    # down to the contig size.
    if fragment_length > len(seq) and ref_circular[contig]:
        if args.small_plasmid_bias:
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
    repeat_length = random.randint(1, 5)
    repeat_count = int(round(fragment_length / repeat_length))
    return get_random_sequence(repeat_length) * repeat_count


def sequence_fragment(fragment, target_identity, error_model, qscore_model):

    # Buffer the fragment a bit so errors can be added to the first and last bases.
    k_size = error_model.kmer_size
    fragment = get_random_sequence(k_size) + fragment + get_random_sequence(k_size)
    frag_len = len(fragment)

    # A list to hold the bases for the errors-added fragment. Note that these values can be ''
    # (meaning the base was deleted) or more than one base (meaning there was an insertion).
    new_fragment_bases = [x for x in fragment]

    errors = 0.0
    change_count = 0

    max_kmer_index = len(new_fragment_bases) - 1 - k_size
    while True:
        # If we have changed almost every base in the fragment, then we can give up (the identity
        # is about as low as we can make it). This is likely to only happen when the target
        # identity is very low (below 60%).
        if change_count > 0.9 * frag_len:
            break

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
                    new_errors = 1
                else:  # insertion
                    new_errors = len(new_base) - 1

                # As the identity gets lower, adding errors has less effect (presumably because
                # adding an error can shift the alignment in a way that makes the overall identity
                # no worse or even better). So we scale our new error count down a bit using our
                # current estimate of the identity.
                errors += new_errors * (estimated_identity ** 1.5)

                # Every now and then we actually align a piece of the new sequence to its original
                # to improve our estimate of the read's identity.
                if change_count % settings.ALIGNMENT_INTERVAL == 0:

                    # If the sequence is short enough, we align the whole thing and get an exact
                    # identity.
                    if frag_len <= settings.ALIGNMENT_SIZE:
                        cigar = edlib.align(fragment, ''.join(new_fragment_bases),
                                            task='path')['cigar']
                        actual_identity = identity_from_edlib_cigar(cigar)
                        errors = (1.0 - actual_identity) * frag_len

                    # If the sequence is longer, we align a random part of the sequence and use
                    # the result to update the error estimate.
                    else:
                        pos = random.randint(0, frag_len - settings.ALIGNMENT_SIZE)
                        pos2 = pos+settings.ALIGNMENT_SIZE
                        cigar = edlib.align(fragment[pos:pos2],
                                            ''.join(new_fragment_bases[pos:pos2]),
                                            task='path')['cigar']
                        actual_identity = identity_from_edlib_cigar(cigar)
                        estimated_errors = (1.0 - actual_identity) * frag_len
                        weight = settings.ALIGNMENT_SIZE / frag_len
                        errors = (estimated_errors * weight) + (errors * (1-weight))

    start_trim = len(''.join(new_fragment_bases[:k_size]))
    end_trim = len(''.join(new_fragment_bases[-k_size:]))

    seq = ''.join(new_fragment_bases)
    qual = get_qscores(seq, fragment, qscore_model)
    assert(len(seq) == len(qual))

    seq = seq[start_trim:-end_trim]
    qual = qual[start_trim:-end_trim]

    return seq, qual


def get_qscores(seq, frag, qscore_model):
    assert len(seq) > 0

    # TODO: I fear this full sequence alignment will be slow for long and inaccurate sequences.
    #       Can I break it into chunks for better performance?
    cigar = edlib.align(seq, frag, task='path')['cigar']
    aligned_seq, aligned_frag, full_cigar = align_sequences_from_edlib_cigar(seq, frag, cigar)
    aligned_len = len(aligned_seq)
    qscores = []
    for i in range(aligned_len):
        if full_cigar[i] == 'D':
            continue
        start, end = i, i
        partial_cigar = full_cigar[i]
        k_size = 1
        while start > 0 and end < aligned_len - 1 and k_size < qscore_model.kmer_size:
            start -= 1
            while full_cigar[start] == 'D':
                start -= 1
            end += 1
            while full_cigar[end] == 'D':
                end += 1
            partial_cigar = full_cigar[start:end+1]
            k_size = len(partial_cigar.replace('D', ''))
            assert k_size % 2 == 1  # should be an odd length k-mer
            if k_size >= qscore_model.kmer_size:
                break
        assert k_size <= qscore_model.kmer_size
        q = qscore_model.get_qscore(partial_cigar)
        qscores.append(q)

    return ''.join(qscores)


def get_start_adapter(rate, amount, adapter):
    if not adapter:
        return ''
    if random_chance(rate):
        if amount == 1.0:
            return adapter
        adapter_frag_length = get_adapter_frag_length(amount, adapter)
        start_pos = len(adapter) - adapter_frag_length
        return adapter[start_pos:]
    return ''


def get_end_adapter(rate, amount, adapter):
    if not adapter:
        return ''
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


def print_glitch_summary(glitch_rate, glitch_size, glitch_skip):
    print('', file=sys.stderr)
    if glitch_rate == 0:
        print('Reads will have no glitches', file=sys.stderr)
    else:
        print('Read glitches:', file=sys.stderr)
        print('  rate (mean distance between glitches): {}'.format(float_to_str(glitch_rate)),
              file=sys.stderr)
        print('  size (mean length of random sequence): {}'.format(float_to_str(glitch_size)),
              file=sys.stderr)
        print('  skip (mean sequence lost per glitch):  {}'.format(float_to_str(glitch_skip)),
              file=sys.stderr)


def add_glitches(fragment, glitch_rate, glitch_size, glitch_skip):
    if glitch_rate == 0:
        return fragment
    i = 0
    new_fragment = []
    while True:
        dist_to_glitch = np.random.geometric(p=1/glitch_rate)
        new_fragment.append(fragment[i:i + dist_to_glitch])
        i += dist_to_glitch
        if i >= len(fragment):
            break

        # Add a glitch!
        if glitch_size > 0:
            new_fragment.append(get_random_sequence(np.random.geometric(p=1/glitch_size)))
        if glitch_skip > 0:
            i += np.random.geometric(p=1/glitch_skip)
        if i >= len(fragment):
            break

    return ''.join(new_fragment)


def print_progress(bp):
    print('\rGenerating reads: {:,} bp'.format(bp), file=sys.stderr, flush=True, end='')
