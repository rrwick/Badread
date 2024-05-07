"""
This module contains code for Badread's error_model subcommand, which is used to generate an error
model using a reference, reads and read-to-reference alignments. This is only needed if users want
to make their own error model instead of using one of the pre-built ones that come with Badread.

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

import collections
import edlib
import itertools
import os
import pathlib
import random
import re
import sys
from .alignment import load_alignments, align_sequences
from .misc import load_fasta, load_fastq, reverse_complement, random_chance, get_random_base, \
    get_random_different_base, get_open_func, check_alignment_matches_read_and_refs, only_acgt


def make_error_model(args, output=sys.stderr, dot_interval=1000):
    refs, _, _ = load_fasta(args.reference)
    reads = load_fastq(args.reads, output=output)
    alignments = load_alignments(args.alignment, args.max_alignments, output=output)
    if len(alignments) == 0:
        sys.exit('Error: no usable alignments')

    kmer_list = [''.join(x) for x in itertools.product('ACGT', repeat=args.k_size)]
    kmer_alternatives = {x: collections.defaultdict(int) for x in kmer_list}

    i = 0
    print('Processing alignments', end='', file=output, flush=True)
    for a in alignments:
        check_alignment_matches_read_and_refs(a, reads, refs)
        read_seq, read_qual = (x[a.read_start:a.read_end] for x in reads[a.read_name])
        ref_seq = refs[a.ref_name][a.ref_start:a.ref_end]

        if a.strand == '-':
            ref_seq = reverse_complement(ref_seq)
        aligned_read_seq, _, aligned_ref_seq, _ = align_sequences(read_seq, read_qual, ref_seq, a)
        start, end = 0, 0
        while True:
            if end > len(aligned_ref_seq):
                break
            ref_kmer = aligned_ref_seq[start:end].replace('-', '')
            if len(ref_kmer) < args.k_size:
                end += 1
                continue
            assert len(ref_kmer) == args.k_size
            read_kmer = aligned_read_seq[start:end].replace('-', '')
            if len(read_kmer) > 1 and ref_kmer[0] == read_kmer[0] and \
                    ref_kmer[-1] == read_kmer[-1] and only_acgt(ref_kmer) and only_acgt(read_kmer):
                kmer_alternatives[ref_kmer][read_kmer] += 1
            start += 1
            while aligned_ref_seq[start] == '-':
                start += 1
            end += 1
        i += 1
        if i % dot_interval == 0:
            print('.', end='', file=output, flush=True)
    print('', file=output, flush=True)

    for kmer in kmer_list:
        alternatives = kmer_alternatives[kmer]
        if len(alternatives) == 0:
            continue
        total = sum(alternatives.values())
        print(f'{kmer},{alternatives[kmer] / total:.6f}', end=';')
        alt_fracs = [(alt_k, count/total) for alt_k, count in alternatives.items() if alt_k != kmer]
        alt_fracs = sorted(alt_fracs, reverse=True, key=lambda x: x[1])
        for k, frac in alt_fracs[:args.max_alt]:
            print(f'{k},{frac:.6f}', end=';')
        print()


class ErrorModel(object):

    def __init__(self, model_type_or_filename, output=sys.stderr):
        self.kmer_size = None
        self.alternatives = {}
        self.probabilities = {}
        this_script_dir = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))

        if model_type_or_filename == 'random':
            print('\nUsing a random error model', file=output)
            self.type = 'random'
            self.kmer_size = 1
        elif model_type_or_filename == 'nanopore2018':
            self.load_from_file(str(this_script_dir / 'error_models' / 'nanopore2018.gz'), output)
        elif model_type_or_filename == 'nanopore2020':
            self.load_from_file(str(this_script_dir / 'error_models' / 'nanopore2020.gz'), output)
        elif model_type_or_filename == 'nanopore2023':
            self.load_from_file(str(this_script_dir / 'error_models' / 'nanopore2023.gz'), output)
        elif model_type_or_filename == 'pacbio2016':
            self.load_from_file(str(this_script_dir / 'error_models' / 'pacbio2016.gz'), output)
        elif model_type_or_filename == 'pacbio2021':
            self.load_from_file(str(this_script_dir / 'error_models' / 'pacbio2021.gz'), output)
        else:
            self.load_from_file(model_type_or_filename, output)

    def load_from_file(self, filename, output):
        print('\nLoading error model from {}'.format(filename), file=output)
        self.type = 'model'
        count = 0
        with get_open_func(filename)(filename, 'rt') as model_file:
            for line in model_file:
                kmer = line.split(',', 1)[0]
                print('\r  ' + kmer, file=output, end='')

                # All k-mers in the model must be the same size.
                if self.kmer_size is None:
                    self.kmer_size = len(kmer)
                else:
                    assert self.kmer_size == len(kmer)

                alternatives = [x.split(',') for x in line.strip().split(';') if x]
                assert alternatives[0][0] == kmer

                self.alternatives[kmer] = [align_kmers(kmer, x[0]) for x in alternatives]
                self.probabilities[kmer] = [float(x[1]) for x in alternatives]
                count += 1
        print(f'\r  done: loaded error distributions for {count} {self.kmer_size}-mers',
              file=output)

    def add_errors_to_kmer(self, kmer):
        """
        Takes a k-mer and returns a (possibly) mutated version of the k-mer, along with the edit
        distance.
        """
        if self.type == 'random':
            return add_one_random_change(kmer)

        if kmer not in self.alternatives:
            return add_one_random_change(kmer)

        alts = self.alternatives[kmer]
        probs = self.probabilities[kmer]

        # The model probabilities for alternate k-mers should total to 1 or a bit less than 1.
        # If less, then the remaining probability is given to random change.
        random_change_prob = 1.0 - sum(probs)
        if random_change_prob > 0.0:
            alts.append(None)
            probs.append(random_change_prob)

        alt = random.choices(alts, weights=probs)[0]
        if alt is None:
            return add_one_random_change(kmer)
        else:
            return alt


def add_one_random_change(kmer):
    result = [x for x in kmer]  # Change 'ACGT' to ['A', 'C', 'G', 'T']
    error_type = random.choice(['s', 'i', 'd'])
    error_pos = random.randint(0, len(kmer) - 1)
    if error_type == 's':  # substitution
        result[error_pos] = get_random_different_base(result[error_pos])
    elif error_type == 'i':  # insertion
        if random_chance(0.5):
            result[error_pos] = result[error_pos] + get_random_base()
        else:
            result[error_pos] = get_random_base() + result[error_pos]
    else:  # deletion
        result[error_pos] = ''
    return result


def align_kmers(kmer, alt):
    """
    This function aligns a k-mer to its alternative, returning a list of strings that join to make
    the alternative, but positioned according to the k-mer.
    Examples:
        align_kmers('ACGT', 'ACGTT') -> ['A', 'C', 'GT', 'T']
        align_kmers('ACGT', 'ACT') -> ['A', 'C', '', 'T']
    """
    # The reference k-mer must be at least a 3-mer (but the alt can be a 2-mer).
    assert len(kmer) > 2
    assert len(alt) > 1

    result = [kmer[0]] + [None] * (len(kmer) - 2) + [kmer[-1]]

    # The alternative k-mer should match at the start and end, so we can trim the sequences down by
    # one before aligning them.
    assert kmer[0] == alt[0] and kmer[-1] == alt[-1]
    kmer, alt = kmer[1:-1], alt[1:-1]

    # Edlib doesn't seem to like it when the first sequence is empty.
    if len(alt) == 0:
        cigar = '{}D'.format(len(kmer))
    else:
        cigar = edlib.align(alt, kmer, task='path')['cigar']

    cigar_parts = re.findall(r'\d+[IDX=]', cigar)
    kmer_pos, alt_pos = 0, 0
    for c in cigar_parts:
        cigar_type = c[-1]
        cigar_size = int(c[:-1])
        if cigar_type == '=' or cigar_type == 'X':
            for i in range(cigar_size):
                result[kmer_pos+1] = alt[alt_pos]
                alt_pos += 1
                kmer_pos += 1
        elif cigar_type == 'D':
            for i in range(cigar_size):
                result[kmer_pos+1] = ''
                kmer_pos += 1
        else:
            assert cigar_type == 'I'
            result[kmer_pos] += alt[alt_pos:alt_pos+cigar_size]
            alt_pos += cigar_size

    # If the insertion landed on the first base, shift it over to the second base (ensures that
    # the first and last base are always the same).
    if len(result[0]) == 2:
        first_base, inserted_base = result[0]
        result[0] = first_base
        result[1] = inserted_base + result[1]
    return result
