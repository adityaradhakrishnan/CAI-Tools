# Given a UCSC Genome Browser Table Output and a genome
# generate the sequence of all listed transcripts.

from string import maketrans
import math
import numpy as np

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
			self.expression = 0.0
			
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
		
		self.length    = len(Sequence)	
		self.sequence  = Sequence
		CodonLists     = []
		CodonLists.append(map(''.join, zip(*[iter(Sequence)]*3)))
		CodonLists.append(map(''.join, zip(*[iter(Sequence[1:-2])]*3)))
		CodonLists.append(map(''.join, zip(*[iter(Sequence[2:-1])]*3)))
		self.codonlist = CodonLists
		
	def ProfilingDensity(self, DensityDict, Samples):
		DensityWrapper = []
		RPKMWrapper    = []
		
		for Idx in Samples:
			Density = np.array([])	
			for IdxN in xrange(0,len(self.exstart)):
				if self.strand == '+':
					Density = np.concatenate((Density, DensityDict[Idx][0][self.chr][self.exstart[IdxN]:self.exstop[IdxN]]))
				else:
					Density = np.concatenate((Density, DensityDict[Idx][1][self.chr][self.exstart[len(self.exstart) - IdxN - 1]:self.exstop[len(self.exstart) - IdxN - 1]][::-1]))
					
			DensityWrapper.append(Density)
			RPKMWrapper.append(np.sum(Density)/self.length*1000)
			
		self.profdensity = DensityWrapper
		self.profRPKM    = RPKMWrapper	
	
	def RNASeqDensity(self, DensityDict, Samples):
		DensityWrapper = []
		RPKMWrapper    = []

		for Idx in Samples:
			Density = np.array([])	
			for IdxN in xrange(0,len(self.exstart)):
				if self.strand == '+':
					Density = np.concatenate((Density, DensityDict[Idx][0][self.chr][self.exstart[IdxN]:self.exstop[IdxN]]))
				else:
					Density = np.concatenate((Density, DensityDict[Idx][1][self.chr][self.exstart[len(self.exstart) - IdxN - 1]:self.exstop[len(self.exstart) - IdxN - 1]][::-1]))

			DensityWrapper.append(Density)
			RPKMWrapper.append(np.sum(Density)/self.length*1000)

		self.mRNAdensity = DensityWrapper
		self.mRNARPKM    = RPKMWrapper		
				
	def CalculateCAI(self,CAIDict):
		CAIList = []
		
		for IdxF in self.codonlist:
			CAI      = 0.0
			for Idx in IdxF:
				CAI += math.log(CAIDict[Idx],2)
				
			CAIList.append(2**(CAI/len(IdxF)))
			
		self.CAI    = CAIList
		CAIStrList  = []

		for Idx in self.CAI:
			Val = np.around(Idx,decimals=1)
			if Val < 0.3:
				Val = 0.3
			else:
				Val = float(str(Val)[0:3])
			CAIStrList.append(Val)

		self.CAIStr = CAIStrList

	def CAIWorstRegion(self,CAIDict):
		CAIWorst = 1.0
		PosWorst = 0

		for Idx in xrange(len(self.codonlist[0])-10):
			CAIList = self.codonlist[0][Idx:Idx+10]
			CAI     = 0.0
			for IdxL in CAIList:
				CAI += math.log(CAIDict[IdxL],2)

			Val     = 2**(CAI/10)

			if Val < CAIWorst:
				CAIWorst = Val
				PosWorst = Idx

		self.CAIWorst = CAIWorst
		self.PosWorst = PosWorst

	def CAIBestRegion(self,CAIDict):
		CAIBest = 0.0
		PosBest = 0

		for Idx in xrange(len(self.codonlist[0])-10):
			CAIList = self.codonlist[0][Idx:Idx+10]
			CAI     = 0.0
			for IdxL in CAIList:
				CAI += math.log(CAIDict[IdxL],2)

			Val     = 2**(CAI/10)

			if Val > CAIBest:
				CAIBest = Val
				PosBest = Idx

		self.CAIBest = CAIBest
		self.PosBest = PosBest

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
		
