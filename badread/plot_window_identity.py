"""
This module contains code for Badread's plot subcommand, which is used to visualise read identities
over a sliding window. This functionality is mainly for debugging purposes - most users probably
won't use it.

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

import matplotlib
import matplotlib.pyplot as plt
import sys

from .alignment import load_alignments, align_sequences
from .misc import load_fasta, load_fastq, reverse_complement
from .qscore_model import qscore_char_to_val


def plot_window_identity(args, output=sys.stdout):
    reads = load_fastq(args.reads, output=output)
    refs, _, _ = load_fasta(args.reference)
    alignments = load_alignments(args.alignment, output=output)

    for a in alignments:
        print(a)
        read_seq, read_qual = (x[a.read_start:a.read_end] for x in reads[a.read_name])
        ref_seq = refs[a.ref_name][a.ref_start:a.ref_end]
        if a.strand == '-':
            ref_seq = reverse_complement(ref_seq)
        _, _, _, errors_per_read_pos = align_sequences(read_seq, read_qual, ref_seq, a)
        positions, identities = get_window_means(errors_per_read_pos, args.window, a.read_start,
                                                 convert_to_identity=True)

        if args.qual:
            read_qual = [qscore_char_to_val(q) for q in read_qual]
            _, qualities = get_window_means(read_qual, args.window, a.read_start,
                                            convert_to_identity=False)
        else:
            qualities = None

        if not args.no_plot:
            plot_one_alignment(positions, identities, qualities, args.window, a,
                               len(reads[a.read_name][0]))


def get_window_means(errors_per_read_pos, window_size, read_start, convert_to_identity=True):
    positions, means = [], []
    window_sum = sum(errors_per_read_pos[:window_size])
    for i in range(len(errors_per_read_pos) - window_size):
        window_start = i
        window_end = i + window_size
        window_centre = i + (window_size // 2)

        if convert_to_identity:
            means.append(100.0 * (1.0 - window_sum / window_size))
        else:
            means.append(window_sum / window_size)
        positions.append(read_start + window_centre)

        window_sum -= errors_per_read_pos[window_start]
        window_sum += errors_per_read_pos[window_end]
    return positions, means


class MyAxes(matplotlib.axes.Axes):
    name = 'MyAxes'

    def drag_pan(self, button, _, x, y):
        matplotlib.axes.Axes.drag_pan(self, button, 'x', x, y)  # pretend key=='x'


matplotlib.projections.register_projection(MyAxes)


def plot_one_alignment(positions, identities, qualities, window_size, alignment, read_length):
    fig, ax1 = plt.subplots(1, 1, figsize=(12, 3), subplot_kw={'projection': 'MyAxes'})
    ax1.plot(positions, identities, '-', color='#8F0505')

    plt.ylabel(f'% identity ({window_size} bp windows)')
    plt.title(f'{alignment.read_name} ({read_length} bp, '
              f'{alignment.percent_identity:.1f}% identity)')
    ax1.set_xlim([0, 10000])
    ax1.set_ylim([50, 100])

    if qualities is not None:
        ax2 = ax1.twinx()
        ax2.plot(positions, qualities, '-', color='#05058F')
        ax2.set_ylim([5, 25])

    fig.canvas.manager.toolbar.pan()
    plt.show()
