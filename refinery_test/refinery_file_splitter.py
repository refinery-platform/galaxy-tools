#!/usr/bin/env python

'''
Test tool for splitting output files from Refinery Test Tools

@author: Scott Ouellette

Input: one text file
Output: N output files based on the amount of input files that got concatenated from Refinery test tool runs

Requires Python v2.7

'''

import re
import argparse


def main(args):
    create_many_files(args.input_file)


def create_many_files(input_file):
    # Split file's content when we see data that wasn't added by test tool runs
    file_content = re.split("Output.*|Input.*", input_file.read())

    sanitized_data = [
        data.lstrip("\n") for data in file_content if data.rstrip("\n")]

    # Create N ouput files based on the number of inputs run through test tools
    for num, file_content in enumerate(sanitized_data):
        open("Output file {}.txt".format(num + 1), 'w').write(file_content)


if __name__ == '__main__':
    version = "%(prog)s 0.1"
    description = "Test tool for running workflows on Galaxy platform from Refinery"
    parser = argparse.ArgumentParser(description=description, version=version)

    parser.add_argument('-i', '--in-file', dest='input_file',
                        type=file, metavar='INPUT_FILE', required=True,
                        help='name of the input file')

    # check argument values for errors
    try:
        args = parser.parse_args()
    except IOError as e:
        parser.error(e)

    main(args)
