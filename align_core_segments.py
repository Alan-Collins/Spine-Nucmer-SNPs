#!/usr/bin/env python3

# AUTHOR      :  ALAN COLLINS
# VERSION     :  v0.1
# DATE        :  2021/4/23
# DESCRIPTION :  Perform separate MUSLCE alignments of each corresponding core segment from NUCMER-aligned core genome segments output by spine.pl

import sys
import argparse
import textwrap as _textwrap
from collections import defaultdict



class LineWrapRawTextHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """
    Short function for argparse that wraps text properly when printing to terminal
    """
    def _split_lines(self, text, width):
        text = self._whitespace_matcher.sub(' ', text).strip()
        return _textwrap.wrap(text, width)


def fasta_to_dict(FASTA_file):
    fasta_dict = {}
    with open(FASTA_file, 'r') as f:
        multifasta = f.read()
    f.close()
    fastas = multifasta.split(">")
    trimmed_fastas = []
    for i in fastas:
        if len(i) != 0:
            trimmed_fastas.append(i)

    fastas = trimmed_fastas

    for i in fastas:
        header = i.split("\n")[0]
        seq = "".join(i.split("\n")[1:])
        fasta_dict[header] = seq

    return fasta_dict


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


corresponding_cores_dict = defaultdict(list)

for file in args.delta_files:
    with open(file, 'r') as fin:
        lines = fin.readlines()
        for i in range(2,len(lines),3):
            bits = lines[i].split()
            try:
                backbone = bits[0][1:] # get rid of '>' symbol at index 0
                core_segment = bits[1]
                corresponding_cores_dict[backbone].append(core_segment)
            except:
                print(lines[i])
# print(corresponding_cores_dict)
            




all_cores_dict = {}

for file in args.core_genome_files:
    fasta_dict = fasta_to_dict(file)
    # for fasta header e.g. ERR712459_core_0001_length_819, ERR712459_core is the genome name common to all contigs.
    genome = "_".join(list(fasta_dict.keys())[0].split('_')[:-3]) 
    all_cores_dict[genome] = fasta_dict

