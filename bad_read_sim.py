
import argparse
import numpy as np
import random
import scipy.special
import sys
import uuid


def get_arguments():
    parser = argparse.ArgumentParser(description='Generate long reads with all sorts of problems')

    required_args = parser.add_argument_group('Required arguments',
                                              description='Quantity can be expressed in absolute '
                                                          'or relative values')
    required_args.add_argument('--reference', type=str, required=True,
                               help='Reference FASTA file')
    required_args.add_argument('--quantity', type=str, required=True,
                               help='Either an absolute value (e.g. 250M) or a relative depth '
                                    '(e.g. 25x)')

    length_args = parser.add_argument_group('Read lengths',
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

    id_args = parser.add_argument_group('Read identities',
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

    problem_args = parser.add_argument_group('Read problems',
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

    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    reference = load_reference(args.reference)
    target_size = get_target_size(reference, args.quantity)
    gamma_k, gamma_t = gamma_parameters(args.mean_read_length, args.read_length_stdev)
    beta_a, beta_b = beta_parameters(args.mean_read_identity, args.read_identity_dist_shape,
                                     args.max_read_identity)

    total_size = 0
    while total_size < target_size:

        # Build a fragment out of parts.
        fragment = [get_start_adapter()]
        main_fragment = get_fragment(gamma_k, gamma_t, reference, args.junk_read_rate,
                                     args.random_read_rate)
        fragment.append(main_fragment)

        while random_chance(args.chimera_rate):
            if random_chance(0.25):
                fragment.append(get_end_adapter())
            if random_chance(0.25):
                fragment.append(get_start_adapter())
            fragment.append(get_fragment(gamma_k, gamma_t, reference, args.junk_read_rate,
                                         args.random_read_rate))
        fragment.append(get_end_adapter())
        fragment = ''.join(fragment)

        seq, quals = sequence_fragment(fragment, beta_a, beta_b, args.max_read_identity,
                                       args.glitches, args.skips)
        read_name = uuid.uuid4()

        print('@{}'.format(read_name))
        print(seq)
        print('+')
        print(quals)

        total_size += len(main_fragment)

    # print('', file=sys.stderr)
    # print('N50: {}'.format(get_n50(reads)), file=sys.stderr)
    # print('Mean identity: {:.4f}'.format(get_mean_id(reads)), file=sys.stderr)


def random_chance(chance):
    return random.random() < chance


def load_reference(reference_fasta_filename):
    pass
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO


def get_target_size(reference, quantity):
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO
    return 1000000


def get_fragment(gamma_k, gamma_t, reference, junk_read_rate, random_read_rate):
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO
    fragment_length = int(round(np.random.gamma(gamma_k, gamma_t)))
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO
    return ''


def sequence_fragment(fragment, beta_a, beta_b, beta_max, glitches, skips):
    read_identity = beta_max * np.random.beta(beta_a, beta_b)
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO
    return '', ''


def get_start_adapter():
    # TODO
    # TODO
    # TODO
    # TODO
    return ''


def get_end_adapter():
    # TODO
    # TODO
    # TODO
    # TODO
    return ''


def gamma_parameters(gamma_mean, gamma_stdev):
    # Shape and rate parameterisation:
    gamma_a = (gamma_mean ** 2) / (gamma_stdev ** 2)
    gamma_b = gamma_mean / (gamma_stdev ** 2)

    # Shape and scale parameterisation:
    gamma_k = (gamma_mean ** 2) / (gamma_stdev ** 2)
    gamma_t = (gamma_stdev ** 2) / gamma_mean

    print('', file=sys.stderr)
    print('Read length gamma distribution parameterisation:', file=sys.stderr)
    print('  shape: alpha = k = ' + '%.4e' % gamma_a, file=sys.stderr)
    print('  rate:   beta = ' + '%.4e' % gamma_b, file=sys.stderr)
    print('  scale: theta = ' + '%.4e' % gamma_t, file=sys.stderr)

    print('', file=sys.stderr)
    print('Theoretical read lengths:', file=sys.stderr)
    print('  N10: {} bp'.format(int(round(find_n_value(gamma_a, gamma_b, 10)))), file=sys.stderr)
    print('  N50: {} bp'.format(int(round(find_n_value(gamma_a, gamma_b, 50)))), file=sys.stderr)
    print('  N90: {} bp'.format(int(round(find_n_value(gamma_a, gamma_b, 90)))), file=sys.stderr)

    return gamma_k, gamma_t


def beta_parameters(beta_mean, beta_shape, beta_max):
    beta_a = (beta_mean / beta_max) * (beta_shape**2)
    beta_b = (1 - (beta_mean / beta_max)) * (beta_shape**2)

    print('', file=sys.stderr)
    print('Read identity beta distribution parameterisation:', file=sys.stderr)
    print('  a = ' + '%.4e' % beta_a, file=sys.stderr)
    print('  b = ' + '%.4e' % beta_b, file=sys.stderr)

    return beta_a, beta_b


def find_n_value(a, b, n):
    target = 1.0 - (n / 100.0)
    bottom_range = 0.0
    top_range = 1.0
    while base_distribution_integral(a, b, top_range) < target:
        bottom_range = top_range
        top_range *= 2
    guess = (bottom_range + top_range) / 2.0
    while True:
        integral = base_distribution_integral(a, b, guess)
        if top_range - bottom_range < 0.01:
            return guess
        if guess == target:
            return guess
        elif integral < target:
            bottom_range = guess
            guess = (bottom_range + top_range) / 2.0
        else:  # integral > target:
            top_range = guess
            guess = (bottom_range + top_range) / 2.0


def base_distribution_integral(a, b, x):
    # TODO: this function bombs out if the value of a is too large.
    # Could I use log-gamma functions to avoid this, perhaps?
    g = scipy.special.gamma(a+1)
    h = inc_gamma(a+1, b*x)
    return (g - h) / g


def inc_gamma(a, b):
    """
    SciPy seems to define the incomplete Gamma function a bit differently than WolframAlpha (which
    I used to do the calc), so this function should represent a WolframAlpha incomplete Gamma.
    https://stackoverflow.com/questions/38713199/incomplete-gamma-function-in-scipy
    """
    return scipy.special.gamma(a) * (1-scipy.special.gammainc(a, b))


def get_n50(reads):
    read_lengths = sorted([x[0] for x in reads], reverse=True)
    n50_target = sum(read_lengths) / 2
    n50, sum_so_far = 0, 0
    for read_length in read_lengths:
        sum_so_far += read_length
        if sum_so_far >= n50_target:
            return read_length


def get_mean_id(reads):
    total_length, total_id = 0, 0.0
    for read_length, read_id in reads:
        total_length += read_length
        total_id += (read_length * read_id)
    return total_id / total_length


if __name__ == '__main__':
    main()
