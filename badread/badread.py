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

import argparse
import sys
from .help_formatter import MyParser, MyHelpFormatter
from .version import __version__


def main():
    parser = MyParser(description='Badread',
                      formatter_class=MyHelpFormatter, add_help=False)

    subparsers = parser.add_subparsers(title='Commands', dest='subparser_name')
    simulate_subparser(subparsers)
    model_subparser(subparsers)
    plot_subparser(subparsers)

    longest_choice_name = max(len(c) for c in subparsers.choices)
    subparsers.help = 'R|'
    for choice, choice_parser in subparsers.choices.items():
        padding = ' ' * (longest_choice_name - len(choice))
        subparsers.help += choice + ': ' + padding
        d = choice_parser.description
        subparsers.help += d[0].lower() + d[1:]  # don't capitalise the first letter
        subparsers.help += '\n'

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--version', action='version', version=__version__,
                           help="Show program's version number and exit")

    # If no arguments were used, print the base-level help which lists possible commands.
    if len(sys.argv) == 1:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if args.subparser_name == 'simulate':
        from .simulate import simulate
        simulate(args)

    elif args.subparser_name == 'model':
        from .error_model import make_error_model
        make_error_model(args)

    elif args.subparser_name == 'plot':
        from .plot_window_identity import plot_window_identity
        plot_window_identity(args)


def simulate_subparser(subparsers):
    group = subparsers.add_parser('simulate', description='Simulate bad reads',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('--reference', type=str, required=True,
                               help='Reference FASTA file')
    required_args.add_argument('--quantity', type=str, required=True,
                               help='Either an absolute value (e.g. 250M) or a relative depth '
                                    '(e.g. 25x)')

    length_args = group.add_argument_group('Read lengths',
                                           description='Read fragments are generated with a '
                                                       'gamma distribution:'
                                                       'desmos.com/calculator/rddlqip1if')
    length_args.add_argument('--mean_read_length', type=int, default=10000,
                             help='Mean read length (in bp)')
    length_args.add_argument('--read_length_stdev', type=int, default=9000,
                             help='Read length standard deviation (in bp)')
    length_args.add_argument('--min_read_length', type=int, default=25,
                             help='Regardless of the distribution, no reads shorter than this '
                                  'will be outputted')

    id_args = group.add_argument_group('Read identities',
                                       description='Read identities are generated with a beta '
                                                   'distribution:'
                                                   'https://www.desmos.com/calculator/t03zr2thap')
    id_args.add_argument('--mean_read_identity', type=float, default=0.85,
                         help='Mean read length (in bp)')
    id_args.add_argument('--read_identity_dist_shape', type=int, default=4,
                         help='Shape parameter - large values produce a tighter distribution '
                              'around the mean')
    id_args.add_argument('--max_read_identity', type=float, default=0.95,
                         help='The entire beta distribution is scaled to max out at this value')
    id_args.add_argument('--error_model', type=str,
                         help='If provided, will use this to simulate realistic errors (otherwise '
                              'errors are random)')

    problem_args = group.add_argument_group('Read problems',
                                            description='Ways reads can go wrong')
    problem_args.add_argument('--junk_read_rate', type=float, default=0.02,
                              help='This fraction of reads will be low-complexity junk')
    problem_args.add_argument('--random_read_rate', type=float, default=0.01,
                              help='This fraction of reads will be random sequence')
    problem_args.add_argument('--chimera_rate', type=float, default=0.01,
                              help='Rate at which separate fragments join together')
    problem_args.add_argument('--start_adapters', type=int, default=25,
                              help='Average amount of adapters on starts of reads')
    problem_args.add_argument('--end_adapters', type=int, default=25,
                              help='Average amount of adapters on ends of reads')
    problem_args.add_argument('--glitches', type=str, default='50,8000',
                              help='Read glitch parameters')
    problem_args.add_argument('--skips', type=str, default='10,8000',
                              help='Read skip parameters')


def model_subparser(subparsers):
    group = subparsers.add_parser('model', description='Generate error model',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('--reference', type=str, required=True,
                               help='Reference FASTA file')
    required_args.add_argument('--reads', type=str, required=True,
                               help='FASTQ of real reads')
    required_args.add_argument('--alignment', type=str, required=True,
                               help='PAF alignment of reads aligned to reference')

    required_args = group.add_argument_group('Optional arguments')
    required_args.add_argument('--k_size', type=int, default=6,
                               help='Error model k-mer size')
    required_args.add_argument('--max_alignments', type=int,
                               help='Only use this many alignments when generating error model')
    required_args.add_argument('--max_alt', type=int, default=100,
                               help='Only save up to this many alternatives to each k-mer')


def plot_subparser(subparsers):
    group = subparsers.add_parser('plot', description='Plot read window identities',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('--reference', type=str, required=True,
                               help='Reference FASTA file')
    required_args.add_argument('--reads', type=str, required=True,
                               help='FASTQ of real reads')
    required_args.add_argument('--alignment', type=str, required=True,
                               help='PAF alignment of reads aligned to reference')

    required_args = group.add_argument_group('Optional arguments')
    required_args.add_argument('--window', type=int, default=100,
                               help='Window size in bp')


if __name__ == '__main__':
    main()
