#!/usr/bin/env python3
import matplotlib.pyplot as plt
import sys
import argparse
import textwrap as _textwrap
import pickle
from itertools import combinations
import multiprocessing
import gzip

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

def seq_diff(seqA, seqB):
	diff = 0
	for a,b in zip(seqA, seqB):
		if a != b:
			diff += 1
	return diff

def build_diff_dict(combo):
	diff_dict[combo] = seq_diff(fasta_dict[combo[0]], fasta_dict[combo[1]])

def pool_make_dist_dict(combo_list, threads, chunksize):

	pool = multiprocessing.Pool(processes=int(threads))

	pool.imap_unordered(build_diff_dict, combo_list, chunksize)
	pool.close()
	pool.join()

parser = argparse.ArgumentParser(
	description="Given a SNP matrix file produced by snps2fasta.py script, calculates all pairwise SNP distances between genomes and stores it in a dict object. This object is then saved as a pickle object. Can also plot the distribution of SNP distances as a simple histogram of either the whole range of SNP distances, or just those up to a defined number of SNPs to focus on more related assemblies.",
	formatter_class=LineWrapRawTextHelpFormatter)
parser.add_argument(
	"-f", dest="input_multifasta", required = False,
	help="Specify input aligned multifasta file produced by snps2fasta.py script."
	)
parser.add_argument(
	"-p", dest="input_pickled_dict", required = False,
	help="Specify input pickled SNP distance dict object created by this script."
	)
parser.add_argument(
	"-d", dest="dict_out", required = False, 
	help="(Optional) Specify output filename with .pkl extension for pickled distance dict of format {(Genome1, Genome2) : #SNPs_between_them}"
	)
parser.add_argument(
	"-z", action='store_true',
		help="Optional compression of the output pickled dict or to declare that input pickled dict is gzipped."
	)
parser.add_argument(
	"-t",  dest="threads", type=int, nargs="?", default = 1,
		help="Specify number of threads to use. Default: 1"
	)
parser.add_argument(
	"-g", dest="hist_out", required = False, 
	help="(Optional) Specify filename for histogram of snp distances."
	)
parser.add_argument(
	"-c", dest="low_snp_hist_out", required = False, 
	help="(Optional) Specify filename for histogram of realtively highly related assemblies snp distances."
	)
parser.add_argument(
	"-ct", dest="low_snp_hist_thresh", required = False, type=int,
	help="(Optional) Specify threshold for histogram of realtively highly related assemblies snp distances. Plots a histogram of only the SNP distances below this level"
	)


args = parser.parse_args(sys.argv[1:])

if args.input_multifasta:
	print("reading in fasta")

	fasta_dict = fasta_to_dict(args.input_multifasta)

	diff_dict = multiprocessing.Manager().dict()
	combos = list(combinations(list(fasta_dict.keys()), 2))
	chunksize = len(combos)//args.threads

	pool_make_dist_dict(combos, args.threads, chunksize)

elif args.input_pickled_dict:
	if not args.z:
		print("reading in pickled distace dict")
		with open(args.input_pickled_dict, 'rb') as pklin:
			diff_dict = pickle.load(pklin)
	else:
		print("reading in gzipped pickled distace dict")
		with gzip.open(args.input_pickled_dict, 'rb') as pklin:
			diff_dict = pickle.load(pklin)




if args.dict_out:
	if args.z:
		try:
			with gzip.open(args.dict_out + ".gz", 'wb') as dump_out:
				pickle.dump(dict(diff_dict), dump_out, protocol=pickle.HIGHEST_PROTOCOL)
		except Exception as e:
			print("pickle_save_gzip_error: " + str(e))
	else:
		try:
			with open(args.dict_out, 'wb') as dump_out:
				pickle.dump(dict(diff_dict), dump_out, protocol=pickle.HIGHEST_PROTOCOL)
		except Exception as e:
			print("pickle_save_error: " + str(e))

if args.hist_out:

	plt.rcParams.update({'font.size': 7})
	fig, ax = plt.subplots()
	fig.set_size_inches(10, 5)
	ax.hist(list(diff_dict.values()), density=False, bins=10000)
	plt.xticks(range(0, max(diff_dict.values())+1, 2000), rotation = 90)
	plt.yscale('log', nonpositive='clip')
	plt.savefig(args.hist_out, dpi=300)
	plt.close()

if args.low_snp_hist_out:
	low_snp_list = []
	for v in diff_dict.values():
		if v < args.low_snp_hist_thresh:
			low_snp_list.append(v)

	plt.rcParams.update({'font.size': 7})
	fig, ax = plt.subplots()
	fig.set_size_inches(10, 5)
	ax.hist(low_snp_list, density=False, bins=range(0,args.low_snp_hist_thresh+100,100))
	plt.xticks(range(0, args.low_snp_hist_thresh+1, 100), rotation = 90)
	plt.savefig(args.low_snp_hist_out, dpi=300)

