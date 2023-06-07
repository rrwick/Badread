#!/usr/bin/env python3
"""
This script takes Badread reads and bins them based on which reference sequence they came from.
For example, if Badread was run on a reference containing 5 different sequences, the script will
create five output files, each containing only the reads from one reference sequence. Random,
junk and chimeric reads will not be included in the output files.

Copyright 2023 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Badread

This file is part of Badread. Badread is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Badread is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Badread.
If not, see <https://www.gnu.org/licenses/>.
"""

import argparse
import collections
import gzip
import pathlib
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='Bin Badread reads by reference sequence')
    parser.add_argument('input_fastq', type=pathlib.Path,
                        help='Filename of input Badread FASTQ file')
    parser.add_argument('output_dir', type=pathlib.Path,
                        help='Output directory name')
    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    pathlib.Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    output_files, output_filenames = {}, {}
    input_count, output_counts = 0, collections.defaultdict(int)
    for _, header, sequence, qualities in iterate_fastq(args.input_fastq):
        input_count += 1
        if ' chimera ' in header or ' random_seq ' in header or ' junk_seq ' in header:
            continue
        try:
            ref_name = header.split(' ')[1].split(',')[0]
        except IndexError:
            continue
        if ref_name not in output_files:
            output_filename = args.output_dir / (ref_name + '.fastq')
            output_filenames[ref_name] = output_filename
            output_files[ref_name] = open(output_filename, 'wt')
        output_files[ref_name].write(f'{header}\n{sequence}\n+\n{qualities}\n')
        output_counts[ref_name] += 1
    for f in output_files.values():
        f.close()
    print('\nInput:')
    print(f'  {args.input_fastq}: {input_count} reads\n')
    print('Output:')
    for ref_name, count in output_counts.items():
        print(f'  {output_filenames[ref_name]}: {count} reads')
    print()


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    https://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)
    unknown_file = open(str(filename), 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('\nError: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('\nError: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


def iterate_fastq(filename):
    with get_open_func(filename)(filename, 'rt') as fastq:
        for line in fastq:
            line = line.strip()
            if len(line) == 0:
                continue
            if not line.startswith('@'):
                continue
            name = line[1:].split()[0]
            header = line
            sequence = next(fastq).strip()
            _ = next(fastq)
            qualities = next(fastq).strip()
            yield name, header, sequence, qualities


if __name__ == '__main__':
    main()
