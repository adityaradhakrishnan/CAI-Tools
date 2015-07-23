# Given a UCSC Genome Browser Table Output and a genome
# generate the sequence of all listed transcripts.

from string import maketrans

class Gene(object):
	
	def __init__(self, InString, GODict, Expression, HalfLives):
		Split = InString.split('\t')
		
		self.name     = Split[1]
		self.chr      = Split[2]
		self.strand   = Split[3]
		self.numexon  = int(Split[8])
		self.exstart  = [int(x) for x in Split[9].split(',')[:-1]]
		self.exstop   = [int(x) for x in Split[10].split(',')[:-1]]
		self.exframe  = [int(x) for x in Split[15].split(',')[:-1]]
		
		# Define ontologies
				
		if self.name in GODict:
			self.ontology = GODict[self.name]
		else:
			self.ontology = []
			
		# Define expression
		
		if self.name in Expression:
			self.expression = Expression[self.name]
		else:
			self.expression = 'DNE'
			
		# Define expression

		if self.name in HalfLives:
			self.halflife = HalfLives[self.name]
		else:
			self.halflife = []
			
	def DefineSeq(self,InString):
		Sequence      = ''
		for Idx in xrange(0,len(self.exstart)):
	        	if self.strand == '+':
                		Sequence += InString[self.exstart[Idx]:self.exstop[Idx]]
            		else:
                		Sequence += InString[self.exstart[len(self.exstart) - Idx - 1]:self.exstop[len(self.exstart) - Idx - 1]][::-1].translate(maketrans("TAGC", "ATCG"))
			
		self.sequence  = Sequence
		self.codonlist = map(''.join, zip(*[iter(Sequence)]*3))
		

	def UWBoundaries(self):
		Count = 0
		for Idx in range(0,len(self.codonlist)-1):
			if self.codonlist[Idx][2] == 'T':
				if self.codonlist[Idx+1][0] == 'A':
					Count += 1
				elif self.codonlist[Idx+1][0] == 'T':
					Count += 1
		
		self.uwboundcount = Count
		self.uwboundfrac  = Count*1.0/(len(self.codonlist) - 1)
		
	def UWGenome(self):
		Count = 0
		for Idx in range(0,len(self.sequence)-1):
			if self.sequence[Idx] == 'T':
				if self.sequence[Idx+1] == 'A':
					Count += 1
				elif self.sequence[Idx+1] == 'T':
					Count += 1
		
		self.uwcount = Count
		self.uwfrac  = Count*1.0/(len(self.sequence) - 1)
		
