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
import os
import pathlib
import random
import re
import sys
from .alignment import load_alignments, align_sequences
from .misc import load_fasta, load_fastq, reverse_complement, float_to_str, get_open_func
from . import settings


def get_qscores(seq, frag, qscore_model):
    assert len(seq) > 0

    # TODO: I fear this full sequence alignment will be slow for long and inaccurate sequences.
    #       Can I break it into chunks for better performance?
    cigar = edlib.align(seq, frag, task='path')['cigar']

    aligned_seq, aligned_frag, full_cigar = align_sequences_from_edlib_cigar(seq, frag, cigar)
    unaligned_len = len(seq)
    margins = (qscore_model.kmer_size - 1) // 2
    qscores = []

    seq_pos_to_alignment_pos = {}
    i, j = 0, 0
    for c in full_cigar:
        if c != 'D':
            seq_pos_to_alignment_pos[i] = j
            i += 1
        j += 1

    for i in range(unaligned_len):
        start = i - margins
        end = i + margins
        while start < 0 or end >= unaligned_len:  # pull back to a smaller k-mer near the seq ends
            start += 1
            end -= 1
        start = seq_pos_to_alignment_pos[start]
        end = seq_pos_to_alignment_pos[end]
        partial_cigar = full_cigar[start:end + 1]
        assert not partial_cigar.startswith('D')
        assert not partial_cigar.endswith('D')
        k_size = len(partial_cigar.replace('D', ''))
        assert k_size <= qscore_model.kmer_size
        assert k_size % 2 == 1  # should be an odd length k-mer
        q = qscore_model.get_qscore(partial_cigar)
        qscores.append(q)

    return ''.join(qscores)


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
                qscore = qscore_char_to_val(qscore)

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


class QScoreModel(object):

    def __init__(self, model_type_or_filename, output=sys.stderr):
        self.scores, self.probabilities = {}, {}
        self.kmer_size = 1
        self.type = None
        this_script_dir = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))

        if model_type_or_filename == 'random':
            self.set_up_random_model(output)
        elif model_type_or_filename == 'ideal':
            self.set_up_ideal_model(output)
        elif model_type_or_filename == 'nanopore':
            self.load_from_file(str(this_script_dir / 'qscore_models' / 'nanopore.gz'), output)
        elif model_type_or_filename == 'pacbio':
            self.load_from_file(str(this_script_dir / 'qscore_models' / 'pacbio.gz'), output)
        else:
            self.load_from_file(model_type_or_filename, output)

        # These three cigars must be in the model, as they are the simplest 1-mer cigars and
        # without them the get_qscore method can fail.
        assert '=' in self.scores
        assert 'X' in self.scores
        assert 'I' in self.scores

    def set_up_random_model(self, output):
        print('\nUsing a random qscore model', file=output)
        self.type = 'random'
        for c in ['=', 'X', 'I']:
            self.scores[c], self.probabilities[c] = \
                uniform_dist_scores_and_probs(settings.RANDOM_QSCORE_MIN,
                                              settings.RANDOM_QSCORE_MAX)

    def set_up_ideal_model(self, output):
        print('\nUsing an ideal qscore model', file=output)
        self.type = 'ideal'
        self.scores['==='], self.probabilities['==='] = \
            uniform_dist_scores_and_probs(settings.IDEAL_QSCORE_GOOD_MIN,
                                          settings.IDEAL_QSCORE_GOOD_MAX)
        self.scores['='], self.probabilities['='] = \
            uniform_dist_scores_and_probs(settings.IDEAL_QSCORE_MEDIUM_MIN,
                                          settings.IDEAL_QSCORE_MEDIUM_MAX)
        for c in ['X', 'I']:
            self.scores[c], self.probabilities[c] = \
                uniform_dist_scores_and_probs(settings.IDEAL_QSCORE_BAD_MIN,
                                              settings.IDEAL_QSCORE_BAD_MAX)

    def load_from_file(self, filename, output):
        print('\nLoading qscore model from {}'.format(filename), file=output)
        self.type = 'model'
        last_cigar_len = 0
        with get_open_func(filename)(filename, 'rt') as model_file:
            for line in model_file:
                parts = line.strip().split(';')
                try:
                    if parts[0] == 'overall':
                        continue

                    cigar = parts[0]
                    k = len(cigar.replace('D', ''))
                    if k > self.kmer_size:
                        self.kmer_size = k
                    print('\r  ' + cigar + (' ' * (last_cigar_len - len(cigar))),
                          file=output, end='')
                    last_cigar_len = len(cigar)
                    scores_and_probs = [x.split(':') for x in parts[2].split(',') if x]
                    self.scores[cigar] = [int(x[0]) for x in scores_and_probs]
                    self.probabilities[cigar] = [float(x[1]) for x in scores_and_probs]
                except (IndexError, ValueError):
                    sys.exit('Error: {} does not seem to be a valid qscore model '
                             'file'.format(filename))
            print('\r  done' + ' ' * (last_cigar_len - 4), file=output)

    def get_qscore(self, cigar):
        """
        If the cigar is in the model, then we use it to choose a qscore. If not, then we trim the
        cigar down by 2 (1 off each end) and try again with the simpler cigar.
        """
        while True:
            assert len(cigar.replace('D', '')) % 2 == 1
            if cigar in self.scores:
                scores = self.scores[cigar]
                probs = self.probabilities[cigar]
                qscore = random.choices(scores, weights=probs)[0]
                break
            else:
                cigar = cigar[1:-1].strip('D')
        return qscore_val_to_char(qscore)


def align_sequences_from_edlib_cigar(seq, frag, cigar, gap_char='-'):
    aligned_seq, aligned_frag, full_cigar = [], [], []
    seq_pos, frag_pos = 0, 0
    cigar_parts = re.findall(r'\d+[IDX=]', cigar)
    for c in cigar_parts:
        cigar_type = c[-1]
        cigar_size = int(c[:-1])
        if cigar_type == '=' or cigar_type == 'X':
            aligned_seq.append(seq[seq_pos:seq_pos+cigar_size])
            aligned_frag.append(frag[frag_pos:frag_pos+cigar_size])
            seq_pos += cigar_size
            frag_pos += cigar_size
        elif cigar_type == 'I':
            aligned_seq.append(seq[seq_pos:seq_pos+cigar_size])
            aligned_frag.append(gap_char * cigar_size)
            seq_pos += cigar_size
        elif cigar_type == 'D':
            aligned_seq.append(gap_char * cigar_size)
            aligned_frag.append(frag[frag_pos:frag_pos+cigar_size])
            frag_pos += cigar_size
        full_cigar.append(cigar_type * cigar_size)
    return ''.join(aligned_seq), ''.join(aligned_frag), ''.join(full_cigar)


def uniform_dist_scores_and_probs(bottom_q, top_q):
    count = top_q - bottom_q + 1
    scores = list(range(bottom_q, top_q + 1))
    probabilities = [1 / count] * count
    return scores, probabilities


def qscore_char_to_val(q):
    return ord(q) - 33


def qscore_val_to_char(q):
    return chr(q + 33)