# Given a UCSC Genome Browser Table Output and a genome
# generate the sequence of all listed transcripts.

import pickle
import time

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
			Sequence += InString[self.exstart[Idx]:self.exstop[Idx]]
			
		self.sequence = Sequence

def ImportGenome(InFile):				
	File      = open(InFile)
	ChromDict = {}

	for line in File:
		if line[0] == '>':
			Split = line.split(' ')
			if len(Split) < 7:
				Name = 'chr' + Split[5].split('=')[1][:-2]
				ChromDict[Name] = ''
			else:
				Name = 'chrM'
				ChromDict[Name] = ''
		else:
			ChromDict[Name] += line.rstrip('\n')
	
	File.close()
	
	File = open(InFile.split('.')[0] + '.pckl','w')
	pickle.dump(ChromDict, File)
	File.close()


PickleGenome = open('2-Ref-Genome/Genome.pckl')
Genome       = pickle.load(PickleGenome)
PickleGenome.close()

File   = open('1-CDS-List/Coding.txt')
Genes  = {}

for line in File:
	Name = line.split('\t')[1]
	Genes[Name] = Gene(line)
	Genes[Name].DefineSeq(Genome[Genes[Name].chr])
	
for IdxG in Genes:
	if (Genes[IdxG].exstop[-1] - Genes[IdxG].exstart[0]) < 100:
		print IdxG