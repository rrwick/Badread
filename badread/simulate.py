"""
This module contains code for Badread's simulate subcommand, which is used to generate fake reads.
This is the most important subcommand in Badread, and for many users the only one they'll ever use.

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
    float_to_str, str_is_int, identity_from_edlib_cigar
from .error_model import ErrorModel
from .qscore_model import QScoreModel, get_qscores
from .fragment_lengths import FragmentLengths
from .identities import Identities
from .version import __version__
from . import settings


def simulate(args, output=sys.stderr):
    print_intro(output)
    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)
    ref_seqs, ref_depths, ref_circular = load_reference(args.reference, output)
    rev_comp_ref_seqs = {name: reverse_complement(seq) for name, seq in ref_seqs.items()}
    frag_lengths = FragmentLengths(args.mean_frag_length, args.frag_length_stdev, output)
    adjust_depths(ref_seqs, ref_depths, ref_circular, frag_lengths, args)
    identities = Identities(args.mean_identity, args.identity_stdev, args.max_identity, output)
    error_model = ErrorModel(args.error_model, output)
    qscore_model = QScoreModel(args.qscore_model, output)
    ref_contigs, ref_contig_weights = get_ref_contig_weights(ref_seqs, ref_depths)
    print_glitch_summary(args.glitch_rate, args.glitch_size, args.glitch_skip, output)

    start_adapt_rate, start_adapt_amount = adapter_parameters(args.start_adapter)
    end_adapt_rate, end_adapt_amount = adapter_parameters(args.end_adapter)
    random_start, random_end = build_random_adapters(args)
    print_adapter_summary(start_adapt_rate, start_adapt_amount, args.start_adapter_seq,
                          end_adapt_rate, end_adapt_amount, args.end_adapter_seq,
                          random_start, random_end, output)

    print_other_problem_summary(args, output)
    ref_size = sum(len(x) for x in ref_seqs.values())
    target_size = get_target_size(ref_size, args.quantity)
    print('', file=output)
    print(f'Target read set size: {target_size:,} bp', file=output)

    print('', file=output)
    count, total_size = 0, 0
    print_progress(count, total_size, target_size, output)
    while total_size < target_size:
        fragment, info = build_fragment(frag_lengths, ref_seqs, rev_comp_ref_seqs, ref_contigs,
                                        ref_contig_weights, ref_circular, args, start_adapt_rate,
                                        start_adapt_amount, end_adapt_rate, end_adapt_amount)
        target_identity = identities.get_identity()
        seq, quals, actual_identity, identity_by_qscores = \
            sequence_fragment(fragment, target_identity, error_model, qscore_model)
        if len(seq) == 0:
            continue

        info.append(f'length={len(seq)}')
        info.append(f'error-free_length={len(fragment)}')
        info.append(f'read_identity={actual_identity * 100.0:.3f}%')

        read_name = uuid.UUID(int=random.getrandbits(128))
        info = ' '.join(info)
        print(f'@{read_name} {info}')
        print(seq)
        print('+')
        print(quals)

        total_size += len(seq)
        count += 1
        print_progress(count, total_size, target_size, output)

    print('\n', file=output)


def build_fragment(frag_lengths, ref_seqs, rev_comp_ref_seqs, ref_contigs, ref_contig_weights,
                   ref_circular, args, start_adapt_rate, start_adapt_amount, end_adapt_rate,
                   end_adapt_amount):
    fragment = [get_start_adapter(start_adapt_rate, start_adapt_amount, args.start_adapter_seq)]
    info = []
    frag_seq, frag_info = get_fragment(frag_lengths, ref_seqs, rev_comp_ref_seqs,
                                       ref_contigs, ref_contig_weights, ref_circular, args)
    fragment.append(frag_seq)
    info.append(','.join(frag_info))

    while random_chance(args.chimeras / 100):  # percentage to fraction
        info.append('chimera')
        if random_chance(settings.CHIMERA_END_ADAPTER_CHANCE):
            fragment.append(args.end_adapter_seq)
        if random_chance(settings.CHIMERA_START_ADAPTER_CHANCE):
            fragment.append(args.start_adapter_seq)
        frag_seq, frag_info = get_fragment(frag_lengths, ref_seqs, rev_comp_ref_seqs,
                                           ref_contigs, ref_contig_weights, ref_circular, args)
        fragment.append(frag_seq)
        info.append(','.join(frag_info))
    fragment.append(get_end_adapter(end_adapt_rate, end_adapt_amount, args.end_adapter_seq))
    fragment = ''.join(fragment)
    fragment = add_glitches(fragment, args.glitch_rate, args.glitch_size, args.glitch_skip)

    return fragment, info


def get_ref_contig_weights(ref_seqs, ref_depths):
    ref_contigs = [x[0] for x in ref_depths.items()]
    ref_contig_weights = [x[1] * len(ref_seqs[x[0]]) for x in ref_depths.items()]
    return ref_contigs, ref_contig_weights


def get_target_size(ref_size, quantity):
    try:
        return int(quantity)
    except ValueError:
        pass
    quantity = quantity.lower()
    try:
        last_char = quantity[-1]
        value = float(quantity[:-1])
        if last_char == 'x':
            return int(round(value * ref_size))
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
    for _ in range(1000):
        seq, info = get_real_fragment(fragment_length, ref_seqs, rev_comp_ref_seqs, ref_contigs,
                                      ref_contig_weights, ref_circular)
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
                      ref_contig_weights, ref_circular):

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
        return seq, info

    # If the reference contig is circular and the fragment length is too long, then we fail to get
    # the read.
    if fragment_length > len(seq) and ref_circular[contig]:
        return '', ''

    start_pos = random.randint(0, len(seq)-1)
    end_pos = start_pos + fragment_length

    # The ending position might be past the end of the sequence. If the sequence is linear, we
    # fix this now.
    if not ref_circular[contig] and end_pos > len(seq):
        end_pos = len(seq)

    info.append(f'{start_pos}-{end_pos}')

    # For circular contigs, we may have to loop the read around the contig.
    if ref_circular[contig]:
        if end_pos <= len(seq):
            return seq[start_pos:end_pos], info
        else:
            looped_end_pos = end_pos - len(seq)
            assert looped_end_pos > 0
        return seq[start_pos:] + seq[:looped_end_pos], info

    # For linear contigs, no looping is needed.
    else:
        return seq[start_pos:end_pos], info


def get_junk_fragment(fragment_length):
    repeat_length = random.randint(1, 5)
    repeat_count = int(round(fragment_length / repeat_length)) + 1
    junk_frag = get_random_sequence(repeat_length) * repeat_count
    return junk_frag[:fragment_length]


def sequence_fragment(fragment, target_identity, error_model, qscore_model):

    # Buffer the fragment a bit so errors can be added to the first and last bases.
    k_size = error_model.kmer_size
    fragment = get_random_sequence(k_size) + fragment + get_random_sequence(k_size)
    frag_len = len(fragment)

    # A list to hold the bases for the errors-added fragment. Note that these values can be ''
    # (meaning the base was deleted) or more than one base (meaning there was an insertion).
    new_fragment_bases = [x for x in fragment]

    errors = 0.0
    change_count, loop_count = 0, 0
    max_kmer_index = len(new_fragment_bases) - 1 - k_size
    estimated_errors_needed = frag_len * (1.0 - target_identity)

    while True:
        # If we need less than 1 error (rounded), we can stop immediately.
        if estimated_errors_needed < 0.5:
            break

        # A precaution to make sure we don't get caught in an infinite loop.
        loop_count += 1
        if loop_count > 100 * frag_len:
            break

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
    qual, actual_identity, identity_by_qscores = get_qscores(seq, fragment, qscore_model)
    assert(len(seq) == len(qual))

    seq = seq[start_trim:-end_trim]
    qual = qual[start_trim:-end_trim]

    return seq, qual, actual_identity, identity_by_qscores


def get_start_adapter(rate, amount, adapter):
    if not adapter or rate == 0.0 or amount == 0.0:
        return ''
    if random_chance(rate):
        if amount == 1.0:
            return adapter
        adapter_frag_length = get_adapter_frag_length(amount, adapter)
        start_pos = len(adapter) - adapter_frag_length
        return adapter[start_pos:]
    return ''


def get_end_adapter(rate, amount, adapter):
    if not adapter or rate == 0.0 or amount == 0.0:
        return ''
    if random_chance(rate):
        if amount == 1.0:
            return adapter
        adapter_frag_length = get_adapter_frag_length(amount, adapter)
        return adapter[:adapter_frag_length]
    return ''


def get_adapter_frag_length(amount, adapter):
    beta_a = 2.0 * amount
    beta_b = 2.0 - beta_a
    return round(int(len(adapter) * np.random.beta(beta_a, beta_b)))


def print_glitch_summary(glitch_rate, glitch_size, glitch_skip, output):
    print('', file=output)
    if glitch_rate == 0:
        print('Reads will have no glitches', file=output)
    else:
        print('Read glitches:', file=output)
        print(f'  rate (mean distance between glitches) = {float_to_str(glitch_rate):>5}',
              file=output)
        print(f'  size (mean length of random sequence) = {float_to_str(glitch_size):>5}',
              file=output)
        print(f'  skip (mean sequence lost per glitch)  = {float_to_str(glitch_skip):>5}',
              file=output)


def print_other_problem_summary(args, output):
    print('', file=output)
    print('Other problems:', file=output)
    print(f'  chimera join rate: {args.chimeras}%', file=output)
    print(f'  junk read rate:    {args.junk_reads}%', file=output)
    print(f'  random read rate:  {args.random_reads}%', file=output)


def adapter_parameters(param_str):
    parts = param_str.split(',')
    if len(parts) == 2:
        try:
            return [float(x) / 100 for x in parts]
        except ValueError:
            pass
    sys.exit('Error: adapter parameters must be two comma-separated values between 0 and 1')


def build_random_adapters(args):
    random_start, random_end = False, False
    if str_is_int(args.start_adapter_seq):
        start_len = int(args.start_adapter_seq)
        args.start_adapter_seq = get_random_sequence(start_len)
        random_start = True
    if str_is_int(args.end_adapter_seq):
        end_len = int(args.end_adapter_seq)
        args.end_adapter_seq = get_random_sequence(end_len)
        random_end = True
    return random_start, random_end


def print_adapter_summary(start_rate, start_amount, start_seq, end_rate, end_amount, end_seq,
                          random_start, random_end, output):
    print('', file=output)
    using_start_adapters = (start_seq and start_rate > 0.0 and start_amount > 0.0)
    using_end_adapters = (end_seq and end_rate > 0.0 and end_amount > 0.0)
    if using_start_adapters:
        random_msg = ' (randomly generated)' if random_start else ''
        print('Start adapter:', file=output)
        print(f'  seq: {start_seq}{random_msg}', file=output)
        print(f'  rate:   {start_rate * 100.0:.1f}%', file=output)
        print(f'  amount: {start_amount * 100.0:.1f}%', file=output)
    else:
        print('Start adapter: none', file=output)
    print('', file=output)
    if using_end_adapters:
        random_msg = ' (randomly generated)' if random_end else ''
        print('End adapter:', file=output)
        print(f'  seq: {end_seq}{random_msg}', file=output)
        print(f'  rate:   {end_rate * 100.0:.1f}%', file=output)
        print(f'  amount: {end_amount * 100.0:.1f}%', file=output)
    else:
        print('End adapter: none', file=output)


def add_glitches(fragment, glitch_rate, glitch_size, glitch_skip):
    if glitch_rate == 0:
        return fragment
    i = 0
    new_fragment = []
    while True:
        p = 1 / glitch_rate if glitch_rate > 1 else 1
        dist_to_glitch = np.random.geometric(p=p)
        new_fragment.append(fragment[i:i + dist_to_glitch])
        i += dist_to_glitch
        if i >= len(fragment):
            break

        # Add a glitch!
        if glitch_size > 0:
            p = 1 / glitch_size if glitch_size > 1 else 1
            new_fragment.append(get_random_sequence(np.random.geometric(p=p)))
        if glitch_skip > 0:
            p = 1 / glitch_skip if glitch_skip > 1 else 1
            i += np.random.geometric(p=p)
        if i >= len(fragment):
            break

    return ''.join(new_fragment)


def print_progress(count, bp, target, output):
    plural = ' ' if count == 1 else 's'
    percent = int(1000.0 * bp / target) / 10
    if percent > 100.0:
        percent = 100.0
    print(f'\rSimulating: {count:,} read{plural}  {bp:,} bp  {percent:.1f}%',
          file=output, flush=True, end='')


def load_reference(reference, output):
    print('', file=output)
    print(f'Loading reference from {reference}', file=output)
    ref_seqs, ref_depths, ref_circular = load_fasta(reference)
    plural = '' if len(ref_seqs) == 1 else 's'
    print(f'  {len(ref_seqs):,} contig{plural}:', file=output)
    for contig in ref_seqs:
        circular_linear = 'circular' if ref_circular[contig] else 'linear'
        print(f'    {contig}: {len(ref_seqs[contig]):,} bp, {circular_linear}, '
              f'{ref_depths[contig]:.2f}x depth', file=output)
    if len(ref_seqs) > 1:
        total_size = sum(len(s) for s in ref_seqs.values())
        print(f'  total size: {total_size:,} bp', file=output)
    return ref_seqs, ref_depths, ref_circular


def print_intro(output):
    print('', file=output)
    print(f'Badread v{__version__}', file=output)
    print(f'long read simulation', file=output)


def adjust_depths(ref_seqs, ref_depths, ref_circular, frag_lengths, args):
    sampled_lengths = [frag_lengths.get_fragment_length() for x in range(100000)]
    total = sum(sampled_lengths)
    for ref_name, ref_seq in ref_seqs.items():
        ref_len = len(ref_seq)
        ref_circ = ref_circular[ref_name]

        # Circular plasmids may have to have their depth increased due compensate for misses.
        if not args.small_plasmid_bias and ref_circ:
            passing_total = sum(length for length in sampled_lengths if length <= ref_len)
            if passing_total == 0:
                sys.exit('Error: fragment length distribution incompatible with reference lengths '
                         '- try running with --small_plasmid_bias to avoid this error')
            adjustment = total / passing_total
            ref_depths[ref_name] *= adjustment

        # Linear plasmids may have to have their depth increased due compensate for truncations.
        if not ref_circ:
            passing_total = sum(min(ref_len, length) for length in sampled_lengths)
            adjustment = total / passing_total
            ref_depths[ref_name] *= adjustment
