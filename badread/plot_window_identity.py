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
import gzip
import matplotlib
import matplotlib.pyplot as plt
import re
import sys


def plot_window_identity(args):
    reads = load_fastq(args.reads)
    refs = load_fasta(args.reference)
    alignments = load_alignments(args.alignment)

    for a in alignments:
        print(a)
        read_seq, read_qual = (x[a.read_start:a.read_end] for x in reads[a.read_name])
        ref_seq = refs[a.ref_name][a.ref_start:a.ref_end]
        if a.strand == '-':
            ref_seq = reverse_complement(ref_seq)
        _, _, errors_per_read_pos = align_sequences(read_seq, ref_seq, a)
        positions, identities = get_window_identity(errors_per_read_pos, args.window, a.read_start)
        plot_one_alignment(positions, identities, args.window, a, len(reads[a.read_name][0]))


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('Error: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('Error: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


REV_COMP_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                 'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
                 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
                 'D': 'H', 'H': 'D', 'N': 'N',
                 'r': 'y', 'y': 'r', 's': 's', 'w': 'w',
                 'k': 'm', 'm': 'k', 'b': 'v', 'v': 'b',
                 'd': 'h', 'h': 'd', 'n': 'n',
                 '.': '.', '-': '-', '?': '?'}


def complement_base(base):
    try:
        return REV_COMP_DICT[base]
    except KeyError:
        return 'N'


def reverse_complement(seq):
    return ''.join([complement_base(x) for x in seq][::-1])


def load_fastq(filename):
    reads = {}
    with get_open_func(filename)(filename, 'rb') as fastq:
        for line in fastq:
            stripped_line = line.strip()
            if len(stripped_line) == 0:
                continue
            if not stripped_line.startswith(b'@'):
                continue
            name = stripped_line[1:].split()[0]
            sequence = next(fastq).strip()
            _ = next(fastq)
            qualities = next(fastq).strip()
            reads[name.decode()] = (sequence.decode(), qualities.decode())
    return reads


def load_fasta(filename):
    fasta_seqs = {}
    with get_open_func(filename)(filename, 'rt') as fasta_file:
        name = ''
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    fasta_seqs[name.split()[0]] = ''.join(sequence)
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line)
        if name:
            fasta_seqs[name.split()[0]] = ''.join(sequence)
    return fasta_seqs


class Alignment(object):

    def __init__(self, paf_line):
        line_parts = paf_line.strip().split('\t')

        self.read_name = line_parts[0]
        self.read_start = int(line_parts[2])
        self.read_end = int(line_parts[3])
        self.strand = line_parts[4]

        self.ref_name = line_parts[5]
        self.ref_length = int(line_parts[6])
        self.ref_start = int(line_parts[7])
        self.ref_end = int(line_parts[8])

        self.matching_bases = int(line_parts[9])
        self.num_bases = int(line_parts[10])
        self.percent_identity = 100.0 * self.matching_bases / self.num_bases

        self.cigar, self.alignment_score = None, None
        for part in line_parts:
            if part.startswith('cg:Z:'):
                self.cigar = part[5:]
            if part.startswith('AS:i:'):
                self.alignment_score = int(part[5:])
        if self.cigar is None:
            sys.exit('Error: no CIGAR string found')
        if self.alignment_score is None:
            sys.exit('Error: no alignment score')

        self.max_indel = 0
        self.cigar_parts = re.findall(r'\d+\w', self.cigar)
        for cigar_part in self.cigar_parts:
            num = int(cigar_part[:-1])
            letter = cigar_part[-1]
            if (letter == 'I' or letter == 'D') and num > self.max_indel:
                self.max_indel = num

        # I want the CIGAR in terms of the read, so I need to flip it if it aligned to the other
        # strand of the reference.
        if self.strand == '-':
            self.cigar_parts = self.cigar_parts[::-1]

    def __repr__(self):
        return self.read_name + ':' + str(self.read_start) + '-' + str(self.read_end) + \
               '(' + self.strand + '),' + \
               self.ref_name + ':' + str(self.ref_start) + '-' + str(self.ref_end) + \
               '(' + ('%.3f' % self.percent_identity) + '%)'


def load_alignments(filename):
    all_alignments = collections.defaultdict(list)
    with get_open_func(filename)(filename, 'rt') as paf_file:
        for line in paf_file:
            a = Alignment(line)
            all_alignments[a.read_name].append(a)
    best_alignments = []
    for read_name, alignments in all_alignments.items():
        best = sorted(alignments, key=lambda x: x.alignment_score)[-1]
        if best.num_bases > 1000 and best.percent_identity > 80.0:
            best_alignments.append(best)
    return best_alignments


def align_sequences(read_seq, ref_seq, alignment):
    read, ref = [], []
    read_pos, ref_pos = 0, 0
    errors_per_read_pos = [0] * len(read_seq)
    alignment.insertions, alignment.deletions, alignment.mismatches = 0, 0, 0
    for c in alignment.cigar_parts:
        cigar_type = c[-1]
        cigar_size = int(c[:-1])
        if cigar_type == 'M':
            read.append(read_seq[read_pos:read_pos+cigar_size])
            ref.append(ref_seq[ref_pos:ref_pos+cigar_size])
            for i in range(cigar_size):
                if read_seq[read_pos+i] != ref_seq[ref_pos+i]:
                    errors_per_read_pos[read_pos+i] += 1
                    alignment.mismatches += 1
            read_pos += cigar_size
            ref_pos += cigar_size
        if cigar_type == 'I':
            read.append(read_seq[read_pos:read_pos+cigar_size])
            ref.append('-' * cigar_size)
            for i in range(cigar_size):
                errors_per_read_pos[read_pos+i] += 1
            alignment.insertions += cigar_size
            read_pos += cigar_size
        if cigar_type == 'D':
            read.append('-' * cigar_size)
            ref.append(ref_seq[ref_pos:ref_pos+cigar_size])
            errors_per_read_pos[read_pos] += cigar_size
            alignment.deletions += cigar_size
            ref_pos += cigar_size
    return ''.join(read), ''.join(ref), errors_per_read_pos


def get_window_identity(errors_per_read_pos, window_size, read_start):
    positions, identities = [], []
    window_sum = sum(errors_per_read_pos[:window_size])
    for i in range(len(errors_per_read_pos) - window_size):
        window_start = i
        window_end = i + window_size
        window_centre = i + (window_size // 2)

        identities.append(100.0 * (1.0 - window_sum / window_size))
        positions.append(read_start + window_centre)

        window_sum -= errors_per_read_pos[window_start]
        window_sum += errors_per_read_pos[window_end]
    return positions, identities


class MyAxes(matplotlib.axes.Axes):
    name = "MyAxes"

    def drag_pan(self, button, key, x, y):
        matplotlib.axes.Axes.drag_pan(self, button, 'x', x, y)  # pretend key=='x'


matplotlib.projections.register_projection(MyAxes)


def plot_one_alignment(positions, identities, window_size, alignment, read_length):
    fig = plt.figure(figsize=(12, 4))
    fig.add_subplot(111, projection='MyAxes')
    plt.plot(positions, identities, '-')
    plt.ylabel('% identity ({} bp windows)'.format(window_size))
    plt.title('{} ({} bp, {:.1f}% identity)'.format(alignment.read_name, read_length,
                                                    alignment.percent_identity))
    plt.gca().set_xlim([0,10000])
    plt.gca().set_ylim([50,100])
    fig.canvas.manager.toolbar.pan()
    plt.show()
