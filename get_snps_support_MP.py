#!/usr/bin/env python3

# Made config with  paste <(ls ../NUCMER/*.snps) <(ls ../SPINE/*.core_coords.txt) <(ls ../*scaffolds_filt.fasta) <(ls ../SAMS/*.sam) <(ls ../SAMS/) | awk '{gsub("../SAMS/","",$5)}1 {gsub("_aligned.sam","",$5)}1' | sed 's/ /\t/g' > config.txt


import argparse, sys, re, time, multiprocessing, math, os, gzip, pickle
import snpclasses

# Takes a config file with 3 columns: Path/to/*core_coords.txt	Path/to/*_scaffolds.fasta	Path/to/*.sam

start_time = time.time()

class config_obj():
	"""class to store information about which backbone corresponds to which region of which scaffold/contig. Takes line from config file"""
	def __init__(self, config_line):
		config = config_line.split("\t")
		self.snps		=	config[0]
		self.coords 	= 	config[1]
		self.scaffold 	= 	config[2]
		self.sam 		= 	config[3]
		self.id			=	config[4].strip()

def rounddown(x):
	return int(math.floor(x / 100.0)) * 100


def fasta_to_dict(FASTA_file):
	fasta_dict = {}
	fastas = FASTA_file.split(">")
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

def rev_comp(base):
	rev_comp_dict = { 'A' : 'T', 'C' : 'G', 'T' : 'A', 'G': 'C', 'a' : 't', 'c' : 'g', 't' : 'a', 'g': 'c'}
	if base not in ['A','T','C','G','a','t','c','g']:
		return base
	else:
		return rev_comp_dict[base]


def cigar_mod_read(seq1, qual1, cigar):
	"""
	According to the M, I, and D components of a cigar string, modifies a seq and quality score strings so that they are in register with refernce sequence 
	returns: modified sequence, modified quality score string
	"""
	cig_sects = re.findall('([0-9]+[A-Z]+)', cigar)

	seq2 =''
	qual2 = ''

	count = 0

	for sect in cig_sects:
		letter = sect[-1]
		sect_len = int(sect[:-1])
		
		if letter == 'M':
			seq2 += seq1[count:count+sect_len] # Add corresponding portion of original seq
			qual2 += qual1[count:count+sect_len]
			count += sect_len

		elif letter == 'D':
			seq2 += '*' * sect_len # Add asterisks for each deleted position relative to the reference sequence
			qual2 += '!' * sect_len # Add ! so that when converted to phred-33 number score will be 0 and won't contribute to consensus score calculation


		elif letter == 'I' or letter == 'S':
			count += sect_len

	return seq2, qual2

def get_read_support(sam_dict, contig, position,snp_rev_comp):
	
	snp_bin = rounddown(position)

	bases, qualities, ids = [], [], []

	for entry in sam_dict[contig][snp_bin].values():
		for i in entry:
			if i.pos <= position and position < i.end:
				try:
					if snp_rev_comp:
						bases.append(rev_comp(i.seq[position - i.pos]))
					else:
						bases.append(i.seq[position - i.pos])
					qualities.append(ord(i.qual[position - i.pos])-33)
					ids.append(i.qname)
				except Exception as e: 
					print("Error processing read %s in contig %s at position %i: " + str(e) %(i.qname, contig, position - i.pos))

	return bases, qualities, ids


def calc_consensus_score(called_base, read_bases, read_qualities, contig, position):
	agree_score = 0
	disagree_score = 0

	for i in range(len(read_bases)):
		base = read_bases[i]
		quality = int(read_qualities[i])
		if base == called_base:
			agree_score += quality
		else:
			disagree_score += quality

	try:
		consensus_score = (agree_score - disagree_score) / agree_score
	except Exception as e:
		print("Exception while calculating consensus score: %s" % e)
		print(contig)
		print(str(position))
		print(called_base)
		print(read_bases)
		print(read_qualities)
		consensus_score = 0

	return str("{:0.2f}".format(consensus_score))


def do_MP_SNP_support(config):
	
	print("Processing %s..." % config.id)

	contig_lookups = {}
	rev_contig_lookups = {} # dict to store contigs with snps so we can only import reads aligning to those contigs.
	snps_list = []
	sam_dict = {}
	relevent_contigs = {}


	print("Reading %s snps file..." % config.id)

	with open(config.snps, 'r') as snps_file:
		snps_lines = snps_file.readlines()[5:]

		if len(snps_lines) != 0:
			try:
				for snps_line in snps_lines:
					try:
						snps_list.append(snpclasses.snp(snps_line))
					except Exception as e:
						print('Exception processing: ' + str(e))
					if snps_line.split()[-1] not in relevent_contigs.keys():
						relevent_contigs[snps_line.split()[-1]] = ''
			except Exception as e:
				print('Exception processing %s: ' + str(e) % config.id)
	snps_file.close()


	print("%i snps found for %s" %(len(snps_lines), config.id))
	
	if len(snps_list) != 0:

		with open(config.coords, 'r') as coords_file:
			for config_line in coords_file.readlines()[1:]:
				try:
					contig_lookups[config_line.split()[4]] = snpclasses.coords_lookup(config_line.split())
				except:
					pass
				try:
					if config_line.split()[-1] in relevent_contigs.keys():
						if config_line.split()[0] in rev_contig_lookups.keys(): # store a list of the start and stop coords of each core contig in the relevant assembly contig
							rev_contig_lookups[config_line.split()[0]].append((config_line.split()[4], int(config_line.split()[2]), int(config_line.split()[3])))
						else:
							rev_contig_lookups[config_line.split()[0]] = [(config_line.split()[4], int(config_line.split()[2]), int(config_line.split()[3]))]
				except Exception as e:
					print("coords_error: " + str(e))
		coords_file.close()

		
		print("Reading %s assembly file..." % config.id)

		with open(config.scaffold, 'r') as scaffold_file:
			fasta_dict = fasta_to_dict(scaffold_file.read())
		scaffold_file.close()

		print("Reading %s sam file..." % config.id)

		if args.pickle and os.path.isfile(config.id + "_sam_dict.pickle.gzip"):
			try:
				with gzip.open(config.id + "_sam_dict.pickle.gzip", 'rb') as picklein:
					sam_dict = pickle.load(picklein)
			except Exception as e:
				print("pickle_save_error: " + str(e))

		else:
			with open(config.sam, 'r') as sam:
				for line in sam.readlines():
					if line[0] != "@":
						entry = snpclasses.SAM_data(line)
						try:
							if entry.rname in rev_contig_lookups.keys():
								for i in rev_contig_lookups[entry.rname]:
									if int(entry.pos) < i[2] and int(entry.end) > i[1]: # For the assembly contig to which this read maps, check if the read maps to a region present in one of the core contigs.
										if any([x in entry.cigar for x in ['D','I','S']]):
											entry.mod_seq, entry.mod_qual = cigar_mod_read(entry.seq, entry.qual, entry.cigar)
											entry.refresh()
										span_locs = [x for x in range(rounddown(entry.pos), entry.end, 100)]
										if entry.rname in sam_dict.keys():
											for loc in span_locs:
												if loc in sam_dict[entry.rname].keys():
													if entry.qname in sam_dict[entry.rname][loc].keys():
														if entry.tlen not in [current.tlen for current in sam_dict[entry.rname][loc][entry.qname]]: # Should be positive for one read and negative for the other if they map to the same contig.
															sam_dict[entry.rname][loc][entry.qname].append(entry)
													else:
														sam_dict[entry.rname][loc][entry.qname] = [entry]
												else:
													sam_dict[entry.rname][loc] = {entry.qname:[entry]}
										else:
											sam_dict[entry.rname] = {loc:{entry.qname:[entry]} for loc in span_locs}
						except Exception as e:
							
							print("sam_error: " + str(e))
			if args.pickle:
				try:	
					with gzip.open(config.id + "_sam_dict.pickle.gzip", 'wb') as pickleout:
						pickle.dump(sam_dict, pickleout, protocol=pickle.HIGHEST_PROTOCOL)
				except Exception as e:
					print("pickle_load_error: " + str(e))


	# 	nreads = 0
	# 	for v in sam_dict.values():
	# 		nreads += len(v)

	# 	print("Found %i reads mapping to core contigs of %s in which SNPs were identified." %(nreads, config.id))

		print("Processing %s reads..." % config.id)

		try:
			print_line_list = []

			for snp in snps_list:
				position = snp.position
				snp_bin = rounddown(int(position))
				variant_call = snp.snp
				contig = snp.contig
				ref_base = snp.ref

				if variant_call != '.' and ref_base != '.':
					snp_loc = contig_lookups[contig].corresponding_pos(position)

					if snp.rev_comp:
						scaffold_nuc = rev_comp(fasta_dict[snp_loc[0]][snp_loc[1]-1])
					else:
						scaffold_nuc = fasta_dict[snp_loc[0]][snp_loc[1]-1]

					print_line = [contig, str(position)]

					for i in snp_loc:
						print_line.append(str(i))

					print_line.append(ref_base)
					print_line.append(variant_call)			
					print_line.append(scaffold_nuc)

					base_calls, base_qualities, read_ids = get_read_support(sam_dict, snp_loc[0], snp_loc[1], snp.rev_comp)

					print_line.append(str(len(base_calls)))
					try:
						print_line.append(str(calc_consensus_score(variant_call, base_calls, base_qualities, contig, position)))
					except Exception as e:
						print("Exception calculating consensus score: " + str(e))
					print_line.append(",".join(base_calls))
					print_line.append(",".join([str(x) for x in base_qualities]))
					print_line.append(",".join(read_ids))

					print_line_list.append(print_line)

		except Exception as e:
			
			print("Exception processing reads: " + str(e))

		print("Writing %s snp support file to %s_snp_support.txt..." % (config.id,config.id))

		with open(args.outdir + "/" + config.id + "_snp_support.txt", 'w+') as outfile:
			outfile.write("Core_contig\tCore_position\tAssembly_contig\tAssembly_position\tReference_base\tVariant_call\tBase_in_assembly\tCoverage\tConsensus_score\tRead_base_calls\tRead_qualities\tRead_identifiers\n")
			outfile.write("\n".join(["\t".join(line) for line in print_line_list]))
		outfile.close()

	else:
		print("No SNPs found for %s. Proceeding to next file." % config.id)


def pool_MP_snp_check(config_list, threads, chunksize):

	pool = multiprocessing.Pool(processes=int(threads))

	pool.imap_unordered(do_MP_SNP_support, config_list, chunksize)
	pool.close()
	pool.join()

class LineWrapRawTextHelpFormatter(argparse.RawDescriptionHelpFormatter):
	"""
	Short function for argparse that wraps text properly when printing to terminal
	"""
	def _split_lines(self, text, width):
		text = self._whitespace_matcher.sub(' ', text).strip()
		return _textwrap.wrap(text, width)


parser = argparse.ArgumentParser(
	description="Given a SNP matrix, reads .fastq files, and core_coords.txt, looks up the read support for each SNP in each assembly and reports that information, highlighting any troubling low confidence base-calls.\n\nUsage example: $ get_snps_support.py -config config.txt",
	formatter_class=LineWrapRawTextHelpFormatter)
parser.add_argument(
	"-config",  dest="config_file", required = True,
	help="Specify config file. "
	)
parser.add_argument(
	"-outdir",  dest="outdir", required = True,
	help="Specify output directory. "
	)
parser.add_argument(
	"-threads",  dest="threads", type=int, nargs="?", default = 1,
		help="Specify number of threads to use. Default: 1"
	)
parser.add_argument(
	"-chunksize",  dest="chunksize", type=int, nargs="?", default = 10,
		help="Specify the number of processes per batch to be sent to each thread. Default: 10"
	)
parser.add_argument(
	"-pickle", action='store_true',
		help="Optional storage of reads that map to regions with snps as a pickled dict of format \{ contig : { 1kb_bin : [list of reads] } }"
	)

args = parser.parse_args(sys.argv[1:])

config_files = []

with open(args.config_file, 'r') as config_file:
	for config_line in config_file.readlines():
		config_files.append(config_obj(config_line))
config_file.close()

pool_MP_snp_check(config_files, args.threads, args.chunksize)

print("Total run time: %.2f" %(time.time()-start_time))
