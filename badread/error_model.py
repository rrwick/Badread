"""
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Badread/

This file is part of Badread. Badread is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Badread is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Badread.
If not, see <http://www.gnu.org/licenses/>.
"""

import collections
import itertools
import sys
from .alignment import load_alignments, align_sequences
from .misc import load_fasta, load_fastq, reverse_complement


def make_error_model(args):
    refs = load_fasta(args.reference)
    reads = load_fastq(args.reads)
    alignments = load_alignments(args.alignments, args.max_alignments)

    kmer_list = [''.join(x) for x in itertools.product('ACGT', repeat=args.k_size)]
    kmer_alternatives = {x: collections.defaultdict(int) for x in kmer_list}

    i = 0
    print('Processing alignments', end='', file=sys.stderr, flush=True)
    for a in alignments:
        # print(a)
        read_seq, read_qual = (x[a.read_start:a.read_end] for x in reads[a.read_name])
        ref_seq = refs[a.ref_name][a.ref_start:a.ref_end]
        if a.strand == '-':
            ref_seq = reverse_complement(ref_seq)
        aligned_read_seq, aligned_ref_seq, errors_per_read_pos = align_sequences(read_seq, ref_seq,
                                                                                 a)

        # print(aligned_read_seq)
        # print(aligned_ref_seq)

        start, end = 0, 0
        while True:
            if end > len(aligned_read_seq):
                break
            read_kmer = aligned_read_seq[start:end].replace('-', '')
            if len(read_kmer) < args.k_size:
                end += 1
                continue
            assert len(read_kmer) == args.k_size
            ref_kmer = aligned_ref_seq[start:end].replace('-', '')
            kmer_alternatives[read_kmer][ref_kmer] += 1

            start += 1
            while aligned_read_seq[start] == '-':
                start += 1
            end += 1
        i += 1
        if i % 1000 == 0:
            print('.', end='', file=sys.stderr, flush=True)
    print('', file=sys.stderr, flush=True)

    for kmer in kmer_list:
        alternatives = kmer_alternatives[kmer]
        if len(alternatives) == 0:
            alternatives[kmer] = 1
        total = sum(alternatives.values())
        print('{},{:.6f}'.format(kmer, alternatives[kmer] / total), end=';')
        alt_fracs = [(alt_k, count/total) for alt_k, count in alternatives.items() if alt_k != kmer]
        alt_fracs = sorted(alt_fracs, reverse=True, key=lambda x: x[1])
        for k, frac in alt_fracs[:args.max_alt]:
            print('{},{:.6f}'.format(k, frac), end=';')
        print()
