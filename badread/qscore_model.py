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
import edlib
import sys
from .alignment import load_alignments, align_sequences
from .misc import load_fasta, load_fastq, reverse_complement


def make_qscore_model(args, output=sys.stderr):
    refs, _, _ = load_fasta(args.reference)
    reads = load_fastq(args.reads)
    alignments = load_alignments(args.alignment, args.max_alignments)

    # The k-mer size has to be odd, so there is a middle base from which we can get the qscore.
    assert args.k_size % 2 == 1

    overall_qscores = collections.defaultdict(int)
    per_cigar_qscores = collections.defaultdict(lambda: collections.defaultdict(int))

    i = 0
    print('Processing alignments', end='', file=output, flush=True)
    for a in alignments:
        read_seq, read_qual = (x[a.read_start:a.read_end] for x in reads[a.read_name])
        ref_seq = refs[a.ref_name][a.ref_start:a.ref_end]
        if a.strand == '-':
            ref_seq = reverse_complement(ref_seq)
        aligned_read_seq, aligned_read_qual, aligned_ref_seq, _ = \
            align_sequences(read_seq, read_qual, ref_seq, a, gap_char=' ')
        start, end = 0, 0
        while True:
            if end > len(aligned_read_seq):
                break
            read_kmer = aligned_read_seq[start:end].replace(' ', '')
            if len(read_kmer) < args.k_size:
                end += 1
                continue
            read_kmer_qual = aligned_read_qual[start:end].replace(' ', '')
            assert len(read_kmer) == len(read_kmer_qual) == args.k_size
            ref_kmer = aligned_ref_seq[start:end].replace(' ', '')

            cigar = edlib.align(read_kmer, ref_kmer, task='path')['cigar']
            qscore = read_kmer_qual[(args.k_size - 1) // 2]
            qscore = ord(qscore) - 33

            overall_qscores[qscore] += 1
            per_cigar_qscores[cigar][qscore] += 1

            start += 1
            while aligned_read_seq[start] == ' ':
                start += 1
            end += 1

        i += 1
        if i % 1000 == 0:
            print('.', end='', file=output, flush=True)
    print('', file=output, flush=True)

    print('overall;', end='')
    print_qscore_fractions(overall_qscores)

    for cigar in sorted(per_cigar_qscores.keys(), key=lambda x: sum(per_cigar_qscores[x].values()),
                        reverse=True):
        print('{};'.format(cigar), end='')
        print_qscore_fractions(per_cigar_qscores[cigar])


def print_qscore_fractions(qscores):
    total = sum(qscores.values())
    for q in sorted(qscores.keys()):
        frac = qscores[q] / total
        print('{};{}:{:.6f},'.format(total, q, frac), end='')
    print()


# class ErrorModel(object):
#
#     def __init__(self, model_type_or_filename, output=sys.stderr):
#         self.kmer_size = None
#         self.alternatives = {}
#         self.probabilities = {}
#         self.mutations_by_base = {}
#
#         print('', file=output)
#         if model_type_or_filename == 'random':
#             print('Using a random error model', file=output)
#             self.type = 'random'
#             self.kmer_size = 1
#         elif model_type_or_filename == 'perfect':
#             print('Using a perfect error model (reads will have no base-level errors)',
#                   file=output)
#             self.type = 'perfect'
#             self.kmer_size = 1
#         else:
#             print('Loading error model from {}'.format(model_type_or_filename), file=output)
#             self.type = 'model'
#             with open(model_type_or_filename, 'rt') as model_file:
#                 for line in model_file:
#                     kmer = line.split(',', 1)[0]
#                     print('\r  ' + kmer, file=output, end='')
#
#                     # All k-mers in the model must be the same size.
#                     if self.kmer_size is None:
#                         self.kmer_size = len(kmer)
#                     else:
#                         assert self.kmer_size == len(kmer)
#
#                     alternatives = [x.split(',') for x in line.strip().split(';') if x]
#                     assert alternatives[0][0] == kmer
#
#                     self.alternatives[kmer] = [align_kmers(kmer, x[0]) for x in alternatives]
#                     self.probabilities[kmer] = [float(x[1]) for x in alternatives]
#             print('\r  done' + ' ' * (self.kmer_size - 4),  # spaces to cover up last k-mer
#                   file=output)
#
#     def add_errors_to_kmer(self, kmer):
#         """
#         Takes a k-mer and returns a (possibly) mutated version of the k-mer, along with the edit
#         distance.
#         """
#         if self.type == 'random':
#             return add_one_random_change(kmer)
#         elif self.type == 'perfect':
#             return [x for x in kmer]
#
#         if kmer not in self.alternatives:
#             return add_one_random_change(kmer)
#
#         alts = self.alternatives[kmer]
#         probs = self.probabilities[kmer]
#
#         # The model probabilities for alternate k-mers should total to 1 or a bit less than 1.
#         # If less, then the remaining probability is given to random change.
#         random_change_prob = 1.0 - sum(probs)
#         if random_change_prob > 0.0:
#             alts.append(None)
#             probs.append(random_change_prob)
#
#         alt = random.choices(alts, weights=probs)[0]
#         if alt is None:
#             return add_one_random_change(kmer)
#         else:
#             return alt
#
#
# def add_one_random_change(kmer):
#     result = [x for x in kmer]  # Change 'ACGT' to ['A', 'C', 'G', 'T']
#     error_type = random.choice(['s', 'i', 'd'])
#     error_pos = random.randint(0, len(kmer) - 1)
#     if error_type == 's':  # substitution
#         result[error_pos] = get_random_different_base(result[error_pos])
#     elif error_type == 'i':  # insertion
#         if random_chance(0.5):
#             result[error_pos] = result[error_pos] + get_random_base()
#         else:
#             result[error_pos] = get_random_base() + result[error_pos]
#     else:  # deletion
#         result[error_pos] = ''
#     return result
#
#
# def align_kmers(kmer, alt):
#     """
#     This function aligns a k-mer to its alternative, returning a list of strings that join to make
#     the alternative, but positioned according to the k-mer.
#     Examples:
#         align_kmers('ACGT', 'ACGTT') -> ['A', 'C', 'GT', 'T']
#         align_kmers('ACGT', 'ACT') -> ['A', 'C', '', 'T']
#     """
#     # The reference k-mer must be at least a 3-mer (but the alt can be a 2-mer).
#     assert len(kmer) > 2
#     assert len(alt) > 1
#
#     result = [kmer[0]] + [None] * (len(kmer) - 2) + [kmer[-1]]
#
#     # The alternative k-mer should match at the start and end, so we can trim the sequences down
#     # by one before aligning them.
#     assert kmer[0] == alt[0] and kmer[-1] == alt[-1]
#     kmer, alt = kmer[1:-1], alt[1:-1]
#
#     # Edlib doesn't seem to like it when the first sequence is empty.
#     if len(alt) == 0:
#         cigar = '{}D'.format(len(kmer))
#     else:
#         cigar = edlib.align(alt, kmer, task='path')['cigar']
#
#     cigar_parts = re.findall(r'\d+[IDX=]', cigar)
#     kmer_pos, alt_pos = 0, 0
#     for c in cigar_parts:
#         cigar_type = c[-1]
#         cigar_size = int(c[:-1])
#         if cigar_type == '=' or cigar_type == 'X':
#             for i in range(cigar_size):
#                 result[kmer_pos+1] = alt[alt_pos]
#                 alt_pos += 1
#                 kmer_pos += 1
#         elif cigar_type == 'D':
#             for i in range(cigar_size):
#                 result[kmer_pos+1] = ''
#                 kmer_pos += 1
#         else:
#             assert cigar_type == 'I'
#             result[kmer_pos] += alt[alt_pos:alt_pos+cigar_size]
#             alt_pos += cigar_size
#
#     # If the insertion landed on the first base, shift it over to the second base (ensures that
#     # the first and last base are always the same).
#     if len(result[0]) == 2:
#         first_base, inserted_base = result[0]
#         result[0] = first_base
#         result[1] = inserted_base + result[1]
#     return result
#
#
# def identity_from_edlib_cigar(cigar):
#     matches, alignment_length = 0, 0
#     cigar_parts = re.findall(r'\d+[IDX=]', cigar)
#     for c in cigar_parts:
#         cigar_type = c[-1]
#         cigar_size = int(c[:-1])
#         alignment_length += cigar_size
#         if cigar_type == '=':
#             matches += cigar_size
#     try:
#         return matches / alignment_length
#     except ZeroDivisionError:
#         return 0.0
