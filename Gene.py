# Given a UCSC Genome Browser Table Output and a genome
# generate the sequence of all listed transcripts.

import pickle

from string import maketrans

class Gene(object):
	
	def __init__(self, InString):
		Split = InString.split('\t')
		
		self.name     = Split[1]
		self.chr      = Split[2]
		self.strand   = Split[3]
		self.numexon  = int(Split[8])
		self.exstart  = [int(x) for x in Split[9].split(',')[:-1]]
		self.exstop   = [int(x) for x in Split[10].split(',')[:-1]]
		self.exframe  = [int(x) for x in Split[15].split(',')[:-1]]
		
	def DefineSeq(self,InString):
		Sequence      = ''
		for Idx in xrange(0,len(self.exstart)):
	        	if self.strand == '+':
                		Sequence += InString[self.exstart[Idx]:self.exstop[Idx]]
            		else:
                		Sequence += InString[self.exstart[len(self.exstart) - Idx - 1]:self.exstop[len(self.exstart) - Idx - 1]][::-1].translate(maketrans("TAGC", "ATCG"))
			
		self.sequence  = Sequence
		self.codonlist = map(''.join, zip(*[iter(Sequence)]*3))
