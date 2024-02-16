"""
This module contains the entry point into Badread. It's run either with the badread command (if it
has been installed) or via the badread-runner.py script (does not require installation). It parses
the arguments and executes the subcommands (which are coded in separate files).

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

import argparse
import pathlib
import sys
from .help_formatter import MyParser, MyHelpFormatter
from .version import __version__
from .misc import bold, str_is_int, str_is_dna_sequence
from . import settings


def main(output=sys.stderr):
    check_python_version()
    args = parse_args(sys.argv[1:])

    if args.subparser_name == 'simulate':
        check_simulate_args(args)
        from .simulate import simulate
        simulate(args, output=output)

    elif args.subparser_name == 'error_model':
        from .error_model import make_error_model
        make_error_model(args, output=output)

    elif args.subparser_name == 'qscore_model':
        from .qscore_model import make_qscore_model
        make_qscore_model(args, output=output)

    elif args.subparser_name == 'plot':
        from .plot_window_identity import plot_window_identity
        plot_window_identity(args)


def parse_args(args):
    parser = MyParser(description=bold('Badread: a long read simulator that can imitate many'
                                       'types of read problems'),
                      formatter_class=MyHelpFormatter, add_help=False)

    subparsers = parser.add_subparsers(title='Commands', dest='subparser_name')
    simulate_subparser(subparsers)
    error_model_subparser(subparsers)
    qscore_model_subparser(subparsers)
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
    help_args.add_argument('--version', action='version', version='Badread v' + __version__,
                           help="Show program's version number and exit")

    # If no arguments were used, print the base-level help which lists possible commands.
    if len(args) == 0:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    return parser.parse_args(args)


def simulate_subparser(subparsers):
    group = subparsers.add_parser('simulate', description='Generate fake long reads',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('--reference', type=str, required=True,
                               help='Reference FASTA file (can be gzipped)')
    required_args.add_argument('--quantity', type=str, required=True,
                               help='Either an absolute value (e.g. 250M) or a relative depth '
                                    '(e.g. 25x)')

    sim_args = group.add_argument_group('Simulation parameters',
                                        description='Length and identity and error distributions')
    sim_args.add_argument('--length', type=str, default='15000,13000',
                          help='Fragment length distribution (mean and stdev, '
                               'default: DEFAULT)')
    sim_args.add_argument('--identity', type=str, default='95,99,2.5',
                          help='Sequencing identity distribution (mean,max,stdev for beta '
                               'distribution or mean,stdev for normal qscore distribution, '
                               'default: DEFAULT)')
    sim_args.add_argument('--error_model', type=str, default='nanopore2023',
                          help='Can be "nanopore2018", "nanopore2020", "nanopore2023", '
                               '"pacbio2016", "pacbio2021", "random" or a model filename')
    sim_args.add_argument('--qscore_model', type=str, default='nanopore2023',
                          help='Can be "nanopore2018", "nanopore2020", "nanopore2023", '
                               '"pacbio2016", "pacbio2021", "random", "ideal" or a model filename')
    sim_args.add_argument('--seed', type=int,
                          help='Random number generator seed for deterministic output (default: '
                               'different output each time)')

    problem_args = group.add_argument_group('Adapters',
                                            description='Controls adapter sequences on the start '
                                                        'and end of reads')
    problem_args.add_argument('--start_adapter', type=str, default='90,60',
                              help='Adapter parameters for read starts (rate and amount, '
                                   'default: DEFAULT)')
    problem_args.add_argument('--end_adapter', type=str, default='50,20',
                              help='Adapter parameters for read ends (rate and amount, '
                                   'default: DEFAULT)')
    problem_args.add_argument('--start_adapter_seq', type=str,
                              default='AATGTACTTCGTTCAGTTACGTATTGCT',
                              help='Adapter sequence for read starts')
    problem_args.add_argument('--end_adapter_seq', type=str, default='GCAATACGTAACTGAACGAAGT',
                              help='Adapter sequence for read ends')

    problem_args = group.add_argument_group('Problems',
                                            description='Ways reads can go wrong')
    problem_args.add_argument('--junk_reads', type=float, default=1,
                              help='This percentage of reads will be low-complexity junk')
    problem_args.add_argument('--random_reads', type=float, default=1,
                              help='This percentage of reads will be random sequence')
    problem_args.add_argument('--chimeras', type=float, default=1,
                              help='Percentage at which separate fragments join together')
    problem_args.add_argument('--glitches', type=str, default='10000,25,25',
                              help='Read glitch parameters (rate, size and skip, default: DEFAULT)')
    problem_args.add_argument('--small_plasmid_bias', action='store_true',
                              help='If set, then small circular plasmids are lost when the '
                                   'fragment length is too high (default: small plasmids are '
                                   'included regardless of fragment length)')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version', version='Badread v' + __version__,
                            help="Show program's version number and exit")


def error_model_subparser(subparsers):
    group = subparsers.add_parser('error_model', description='Build a Badread error model',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('--reference', type=str, required=True,
                               help='Reference FASTA file')
    required_args.add_argument('--reads', type=str, required=True,
                               help='FASTQ of real reads')
    required_args.add_argument('--alignment', type=str, required=True,
                               help='PAF alignment of reads aligned to reference')

    required_args = group.add_argument_group('Optional arguments')
    required_args.add_argument('--k_size', type=int, default=7,
                               help='Error model k-mer size')
    required_args.add_argument('--max_alignments', type=int,
                               help='Only use this many alignments when generating error model '
                                    '(default: use all alignments)')
    required_args.add_argument('--max_alt', type=int, default=25,
                               help='Only save up to this many alternatives to each k-mer')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version', version='Badread v' + __version__,
                            help="Show program's version number and exit")


def qscore_model_subparser(subparsers):
    group = subparsers.add_parser('qscore_model', description='Build a Badread qscore model',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('--reference', type=str, required=True,
                               help='Reference FASTA file')
    required_args.add_argument('--reads', type=str, required=True,
                               help='FASTQ of real reads')
    required_args.add_argument('--alignment', type=str, required=True,
                               help='PAF alignment of reads aligned to reference')

    required_args = group.add_argument_group('Optional arguments')
    required_args.add_argument('--k_size', type=int, default=9,
                               help='Qscore model k-mer size (must be odd, default: DEFAULT)')
    required_args.add_argument('--max_alignments', type=int,
                               help='Only use this many alignments when generating qscore model '
                                    '(default: use all alignments)')
    required_args.add_argument('--max_del', type=int, default=6,
                               help='Deletion runs longer than this will be collapsed to reduce '
                                    'the number of possible alignments')
    required_args.add_argument('--min_occur', type=int, default=100,
                               help='CIGARs which occur less than this many times will not be '
                                    'included in the model')
    required_args.add_argument('--max_output', type=int, default=10000,
                               help='The outputted model will be limited to this many lines')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version', version='Badread v' + __version__,
                            help="Show program's version number and exit")


def plot_subparser(subparsers):
    group = subparsers.add_parser('plot', description='View read identities over a sliding window',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('--reference', type=str, required=True,
                               help='Reference FASTA file')
    required_args.add_argument('--reads', type=str, required=True,
                               help='FASTQ of real reads')
    required_args.add_argument('--alignment', type=str, required=True,
                               help='PAF alignment of reads aligned to reference')

    optional_args = group.add_argument_group('Optional arguments')
    optional_args.add_argument('--window', type=int, default=100,
                               help='Window size in bp')
    optional_args.add_argument('--qual', action='store_true',
                               help='Include qscores in plot (default: only show identity)')
    optional_args.add_argument('--no_plot', action='store_true',
                               help='Do not display plots (for testing purposes)')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version', version='Badread v' + __version__,
                            help="Show program's version number and exit")


def check_simulate_args(args):
    if not pathlib.Path(args.reference).is_file():
        sys.exit(f'Error: {args.reference} is not a file')

    error_model_names = ['random', 'nanopore2018', 'nanopore2020', 'nanopore2023',
                         'pacbio2016', 'pacbio2021']
    error_model = args.error_model.lower()
    if error_model not in error_model_names and not pathlib.Path(args.error_model).is_file():
        sys.exit(f'Error: {args.error_model} is not a file\n'
                 f'  --error_model must be from {error_model_names} or a filename')

    qscore_model_names = ['random', 'ideal', 'nanopore2018', 'nanopore2020', 'nanopore2023',
                          'pacbio2016', 'pacbio2021']
    qscore_model = args.qscore_model.lower()
    if qscore_model not in qscore_model_names and not pathlib.Path(args.qscore_model).is_file():
        sys.exit(f'Error: {args.qscore_model} is not a file\n'
                 f'  --qscore_model must be from {qscore_model_names} or a filename')

    if args.chimeras > 50:
        sys.exit('Error: --chimeras cannot be greater than 50')
    if args.junk_reads > 100:
        sys.exit('Error: --junk_reads cannot be greater than 100')
    if args.random_reads > 100:
        sys.exit('Error: --random_reads cannot be greater than 100')
    if args.junk_reads + args.random_reads > 100:
        sys.exit('Error: --junk_reads and --random_reads cannot sum to more than 100')

    try:
        length_parameters = [float(x) for x in args.length.split(',')]
        args.mean_frag_length = length_parameters[0]
        args.frag_length_stdev = length_parameters[1]
    except (ValueError, IndexError):
        sys.exit('Error: could not parse --length values')
    if args.mean_frag_length <= settings.MIN_MEAN_READ_LENGTH:
        sys.exit(f'Error: mean read length must be at least {settings.MIN_MEAN_READ_LENGTH}')
    if args.frag_length_stdev < 0:
        sys.exit('Error: read length stdev cannot be negative')

    try:
        identity_parameters = [float(x) for x in args.identity.split(',')]
        if len(identity_parameters) == 2:
            args.mean_identity = identity_parameters[0]
            args.max_identity = None
            args.identity_stdev = identity_parameters[1]
            check_qscore_identities(args)
        elif len(identity_parameters) == 3:
            args.mean_identity = identity_parameters[0]
            args.max_identity = identity_parameters[1]
            args.identity_stdev = identity_parameters[2]
            check_beta_identities(args)
        else:
            sys.exit('Error: could not parse --identity values')
    except (ValueError, IndexError):
        sys.exit('Error: could not parse --identity values')

    try:
        glitch_parameters = [float(x) for x in args.glitches.split(',')]
        args.glitch_rate = glitch_parameters[0]
        args.glitch_size = glitch_parameters[1]
        args.glitch_skip = glitch_parameters[2]
    except (ValueError, IndexError):
        sys.exit('Error: could not parse --glitches values')
    if args.glitch_rate < 0 or args.glitch_size < 0 or args.glitch_skip < 0:
        sys.exit('Error: --glitches must contain non-negative values')

    if args.start_adapter_seq != '':
        if not str_is_int(args.start_adapter_seq):
            args.start_adapter_seq = args.start_adapter_seq.upper()
            if not str_is_dna_sequence(args.start_adapter_seq):
                sys.exit('Error: --start_adapter_seq must be a DNA sequence or a number')
    if args.end_adapter_seq != '':
        if not str_is_int(args.end_adapter_seq):
            args.end_adapter_seq = args.end_adapter_seq.upper()
            if not str_is_dna_sequence(args.end_adapter_seq):
                sys.exit('Error: --end_adapter_seq must be a DNA sequence or a number')


def check_beta_identities(args):
    if args.mean_identity > 100.0:
        sys.exit('Error: mean read identity cannot be more than 100')
    if args.max_identity > 100.0:
        sys.exit('Error: max read identity cannot be more than 100')
    if args.mean_identity <= settings.MIN_MEAN_READ_IDENTITY:
        sys.exit(f'Error: mean read identity must be at least {settings.MIN_MEAN_READ_IDENTITY}')
    if args.max_identity <= settings.MIN_MEAN_READ_IDENTITY:
        sys.exit(f'Error: max read identity must be at least {settings.MIN_MEAN_READ_IDENTITY}')
    if args.mean_identity > args.max_identity:
        sys.exit(f'Error: mean identity ({args.mean_identity}) cannot be larger than max '
                 f'identity ({args.max_identity})')
    if args.identity_stdev < 0.0:
        sys.exit('Error: read identity stdev cannot be negative')


def check_qscore_identities(args):
    if args.mean_identity <= settings.MIN_MEAN_READ_QSCORE:
        sys.exit(f'Error: mean read identity must be at least {settings.MIN_MEAN_READ_QSCORE}')
    if args.identity_stdev < 0.0:
        sys.exit('Error: read qscore stdev cannot be negative')


def check_python_version():
    if sys.version_info.major < 3 or sys.version_info.minor < 6:
        sys.exit('Error: Badread requires Python 3.6 or later')


if __name__ == '__main__':
    main()
