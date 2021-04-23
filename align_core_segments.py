#!/usr/bin/env python3

# AUTHOR      :  ALAN COLLINS
# VERSION     :  v1
# DATE        :  2021/4/23
# DESCRIPTION :  Perform separate MUSLCE alignments of each corresponding core segment from NUCMER-aligned core genome segments output by spine.pl

import sys
import argparse
import textwrap as _textwrap


class LineWrapRawTextHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """
    Short function for argparse that wraps text properly when printing to terminal
    """
    def _split_lines(self, text, width):
        text = self._whitespace_matcher.sub(' ', text).strip()
        return _textwrap.wrap(text, width)

parser = argparse.ArgumentParser(
    description="Perform separate MUSLCE alignments of each corresponding core segment from NUCMER-aligned core genome segments output by spine.pl",
    formatter_class=LineWrapRawTextHelpFormatter)
parser.add_argument(
    "-d", dest="delta_files", required = True, nargs="+",
    help="Provide a list of .delta files"
    )
parser.add_argument(
    "-c", dest="core_genome_files", required = True, nargs="+",
    help="Provide a list of core genome fasta files"
    )
parser.add_argument(
    "-o", dest="output_dir", required = True,
    help="Provide the path to the directory you want the output alignments stored"
    )


args = parser.parse_args(sys.argv[1:])

