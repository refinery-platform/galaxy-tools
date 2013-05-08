#!/usr/bin/env python

'''
Test tool for running Galaxy workflows from Refinery

@author: Ilya Sytchev

Input: one or more text files
Output: concatenated input file(s) plus annotation

Requires Python v2.7

'''


import argparse
import random
import sys
import time


def main(args):
    time.sleep(args.seconds)

    if check_fail(args.p_fail):
        cleanup(args)
        quit("Processing failed by request", args)

    input = read_files(args.input_files)
    try:
        for out_file in args.output_files:
            out_file.write('\n' + "Output file name: " + out_file.name + '\n')
            out_file.write(input)
    except IOError as e:
        cleanup(args)
        parser.error(e)
    else:
        cleanup(args)


def check_fail(p_fail):
    '''Determine success/failure state given a probability of failure

    '''
    random.seed()
    if random.random() < p_fail:
        return True
    else:
        return False


def read_files(file_list):
    '''Read files from disk into a string

    '''
    str = ''
    for in_file in file_list:
        str += "Input file name: " + in_file.name + '\n'
        str += in_file.read()
    return str


def cleanup(args):
    '''Close all open file handles

    '''
    file_list = []
    if args.input_files:
        file_list.extend(args.input_files)
    if args.output_files:
        file_list.extend(args.output_files)
    for fh in file_list:
        try:
            fh.close()
        except AttributeError:
            continue


def quit(message, args):
    '''Exit and optionally write to stdout and/or stderr

    '''
    if args.stdout:
        sys.stdout.write(message + '\n')
    if args.stderr:
        sys.stderr.write(message + '\n')
    sys.exit(args.exit_code)


if __name__ == '__main__':
    version = "%(prog)s 0.1"
    description = "Test tool for running workflows on Galaxy platform from Refinery"
    parser = argparse.ArgumentParser(description=description, version=version)

    parser.add_argument('-i', '--in-file', dest='input_files', nargs='+',
                        type=file, metavar='INPUT_FILE', required=True,
                        help='name of the input file')
    parser.add_argument('-o', '--out-file', dest='output_files', nargs='+',
                        type=argparse.FileType('w'), metavar='OUTPUT_FILE',
                        required=True, help='name of the output file')
    parser.add_argument('-e', '--exit_code', type=int, default=0,
                        help='code to return on exit, default: %(default)s')
    parser.add_argument('--stdout', action='store_true',
                        help='write a message to stdout')
    parser.add_argument('--stderr', action='store_true',
                        help='write a message to stderr')
    parser.add_argument('-p', '--p-fail', type=float, default=0.0,
                        help='probability of execution failure, default: %(default)s')
    parser.add_argument('-s', '--sleep', dest='seconds', type=int, default=0,
                        metavar='SECONDS',
                        help='number of seconds to sleep, default: %(default)s')

    # check argument values for errors
    try:
        args = parser.parse_args()
    except IOError as e:
        parser.error(e)

    if args.exit_code < 0 or args.exit_code > 255:
        cleanup(args)
        parser.error("Exit code value must be between 0 and 255")

    if args.p_fail < 0.0 or args.p_fail > 1.0:
        cleanup(args)
        parser.error("Probability value must be between 0.0 and 1.0 inclusive")

    main(args)
