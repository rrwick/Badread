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

import collections
import re
import sys
from .alignment import load_alignments, align_sequences
from .misc import load_fasta, load_fastq, reverse_complement, float_to_str


def make_qscore_model(args, output=sys.stderr):
    refs, _, _ = load_fasta(args.reference)
    reads = load_fastq(args.reads)
    alignments = load_alignments(args.alignment, args.max_alignments)

    # The k-mer size has to be odd, so there is a middle base from which we can get the qscore.
    assert args.k_size % 2 == 1

    overall_qscores = collections.defaultdict(int)
    per_cigar_qscores = collections.defaultdict(lambda: collections.defaultdict(int))

    p = re.compile('D{' + str(args.max_del) + ',}')
    max_del = 'D' * args.max_del

    i = 0
    print('Processing alignments', end='', file=output, flush=True)
    for a in alignments:
        if a.read_name not in reads:
            sys.exit('\nError: could not find read {}\n'
                     'are you sure your read file and alignment file match?'.format(a.read_name))
        if a.ref_name not in refs:
            sys.exit('\nError: could not find reference {}\nare you sure '
                     'your reference file and alignment file match?'.format(a.ref_name))

        read_seq, read_qual = (x[a.read_start:a.read_end] for x in reads[a.read_name])
        ref_seq = refs[a.ref_name][a.ref_start:a.ref_end]

        if a.strand == '-':
            ref_seq = reverse_complement(ref_seq)
        aligned_read_seq, aligned_read_qual, aligned_ref_seq, _ = \
            align_sequences(read_seq, read_qual, ref_seq, a, gap_char=' ')

        for k_size in range(1, args.k_size+2, 2):  # Do all odd k-mer sizes up to the setting
            start, end = 0, 0
            while True:
                if end > len(aligned_read_seq):
                    break
                read_kmer = aligned_read_seq[start:end]
                if len(read_kmer.replace(' ', '')) < k_size:
                    end += 1
                    continue
                read_kmer_qual = aligned_read_qual[start:end].replace(' ', '')
                assert len(read_kmer.replace(' ', '')) == len(read_kmer_qual) == k_size
                ref_kmer = aligned_ref_seq[start:end]

                cigar = []
                for i, read_base in enumerate(read_kmer):
                    ref_base = ref_kmer[i]
                    assert read_base != ' ' or ref_base != ' '
                    if read_base == ref_base:
                        cigar.append('=')
                    elif read_base == ' ':
                        cigar.append('D')
                    elif ref_base == ' ':
                        cigar.append('I')
                    else:
                        cigar.append('X')
                cigar = ''.join(cigar)
                assert len(cigar.replace('D', '')) == k_size
                cigar = p.sub(max_del, cigar)

                qscore = read_kmer_qual[(k_size - 1) // 2]
                qscore = ord(qscore) - 33

                if k_size == 1:
                    overall_qscores[qscore] += 1
                per_cigar_qscores[cigar][qscore] += 1

                start += 1
                if start >= len(aligned_read_seq):
                    break
                while aligned_read_seq[start] == ' ':
                    start += 1
                end += 1

        i += 1
        if i % 1000 == 0:
            print('.', end='', file=output, flush=True)
    print('', file=output, flush=True)

    print_qscore_fractions('overall', overall_qscores, 0)
    for cigar in sorted(per_cigar_qscores.keys(),
                        key=lambda x: (len(x.replace('D', '')),
                                       1 / sum(per_cigar_qscores[x].values()))):
        print_qscore_fractions(cigar, per_cigar_qscores[cigar], args.min_occur)


def print_qscore_fractions(cigar, qscores, min_occur):
    total = sum(qscores.values())
    if total < min_occur:
        return
    print('{};'.format(cigar), end='')
    print('{};'.format(total), end='')
    for q in sorted(qscores.keys()):
        frac = qscores[q] / total
        frac_str = float_to_str(frac, decimals=6, trim_zeros=True)
        print('{}:{},'.format(q, frac_str), end='')
    print()
