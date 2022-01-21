#!/usr/bin/env python3

class snp():
	"""Class to store coordinates and sequence of snps"""
	def __init__(self, snps_line):
		bits 			= snps_line.split()
		self.contig 	= bits[-1]
		self.position 	= bits[3]
		self.snp 		= bits[2]
		self.ref		= bits[1]

		if bits[-3] == "-1":
			self.rev_comp   = True
		else:
			self.rev_comp	= False



class SAM_data(object):
	"""stores columns of SAM entry as attributes"""
	def __init__(self, object):
		self.qname = object.split('\t')[0]
		self.flag = object.split('\t')[1]
		self.rname = object.split('\t')[2]
		self.pos = int(object.split('\t')[3])
		self.mapq = int(object.split('\t')[4])
		self.cigar = object.split('\t')[5]
		self.rnext = object.split('\t')[6]
		self.pnext = object.split('\t')[7]
		self.tlen = object.split('\t')[8]
		self.seq = object.split('\t')[9]
		self.qual = object.split('\t')[10]
		self.ln = len(self.seq)
		self.end = self.pos + self.ln
		self.mod_seq = ''
		self.mod_qual = ''
		self.refreshed = False

	def refresh(self):
		if not self.refreshed:
			self.seq = self.mod_seq
			self.qual = self.mod_qual
			self.ln = len(self.seq)
			self.end = self.pos + self.ln
			self.refreshed = True

class Read_padding():
	"""Store information about how much to pad reads when printing them in order to maintain register between lines"""
	def __init__(self):
		self.read_name = ''
		self.seq = ''
		self.qual = ''
		self.cigar = ''
		self.pad = 0 # 0: no pad, 1: left pad, 2: right pad, 3: both sides pad
		self.padl = 0
		self.padr = 0
		self.warning = ''
		
	def report_pad_seq(self):
		for i in self.qual:
			if ord(i) < 53:
				self.warning = '!!!WARNING: Low quality base call in sequence!!!'

		pseq = '-' * self.padl + self.seq + '-' * self.padr 
		pqual = '-' * self.padl + self.qual + '-' * self.padr
		return '\t'.join([self.read_name, pseq, pqual, self.cigar, self.warning])


class coords_lookup():
	"""Class to store corresponding core and assembly contigs and the position in the assembly contig at which the core contig begins. Takes a line from a coords file."""
	def __init__(self, coords):
		self.assembly_contig 	= coords[0]
		self.core_contig 		= coords[4]
		self.start_pos			= int(coords[2])

	def corresponding_pos(self, position):
		"""Given an integer position within the core cotig, return the corresponding position in the assembly contig. Returns both assembly contig ID and position"""
		return (self.assembly_contig, int(position) + self.start_pos - 1) # -1 because it is adjusting the 1-base coordinates to a 0-base system