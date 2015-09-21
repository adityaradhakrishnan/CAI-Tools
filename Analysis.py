from Gene import Gene
import SequenceTools as st
import numpy as np
import pickle
import time

Start = time.time()

def LoadDict(inString, LengthDict):
	Density = {0:{},1:{}}
	
	String = ['-P.wig','-M.wig']
	
	for Idx in range(0,2):
		File = open(inString + String[Idx])
	
		for line in File:
			if line[0] == 'f':
				Chr               = line.strip('\n').split('  ')[1][6:]
				Density[Idx][Chr] = np.array([0.0]*LengthDict[Chr])
				IdxD              = 0
			elif line[0] == 't':
				pass
			else:
				Density[Idx][Chr][IdxD]  = float(line)
				IdxD                    += 1
				
	return Density

# Load the reference genome

PickleNames    = open('1-CDS-List/NameConvert.pckl')
PickleGenome   = open('2-Ref-Genome/Genome.pckl')
PickleGODict   = open('4-Ontologies/GOSlim.pckl')
PickleExpress  = open('5-Expression/Expression.pckl')
PickleHL       = open('6-RNAHalfLives/HalfLives.pckl')
PickleCAICodon = open('7-CAI/Codon-CAI.pckl')


Names         = pickle.load(PickleNames)
Genome        = pickle.load(PickleGenome)
GODict        = pickle.load(PickleGODict)
Expression    = pickle.load(PickleExpress)
HalfLives     = pickle.load(PickleHL)
CAIDict       = pickle.load(PickleCAICodon)
 
PickleNames.close()
PickleGenome.close()
PickleGODict.close()
PickleExpress.close()
PickleHL.close()
PickleCAICodon.close()

# Load in the necessary WIG files:

LengthDict = {
      'chrI':230218,   'chrII':813184,  'chrIII':316620, 'chrIV':1531933,   'chrV':576874,  
	 'chrVI':270161, 'chrVII':1090940, 'chrVIII':562643,  'chrIX':439888,   'chrX':745751, 
     'chrXI':666816, 'chrXII':1078177, 'chrXIII':924431, 'chrXIV':784333, 'chrXV':1091291,
    'chrXVI':948066,     'chrM':85779
}


File  = open('1-CDS-List/Coding.txt')
Genes = {}
    
for line in File:
	Name = line.split('\t')[1]
	Genes[Name] = Gene(line, GODict, Expression, HalfLives)
	Genes[Name].DefineSeq(Genome[Genes[Name].chr])
	Genes[Name].CalculateCAI(CAIDict)
	Genes[Name].CAIWorstRegion(CAIDict)
	Genes[Name].CAIBestRegion(CAIDict)

RNAString  = '../WIGs/RNASeq/Unmod/'
ProfString = '../WIGs/Profiling/'

RNAType    = ['HBGN1ADXX-2-ACTTGA-WT','HBGN1ADXX-1-TGACCA-DhKO','HBGN1ADXX-1-CAGATC-DhOE']
ProfType   = ['FC278-L6-P1-ACTTGA-WT','FC278-L5-P1-TGACCA-DhKO','FC278-L5-P1-CAGATC-DhOE']

#RNASeqSamples = [RNAString + Idx for Idx in RNAType]
ProfSamples   = [ProfString + Idx for Idx in ProfType] 

#RNASeqDensityDict = {}   
ProfDensityDict   = {}   

#for Idx in RNASeqSamples:
	#print Idx
	#RNASeqDensityDict[Idx] = LoadDict(Idx, LengthDict)

for Idx in ProfSamples:
	print Idx
	ProfDensityDict[Idx]   = LoadDict(Idx, LengthDict)	

BestWT  = np.array([0]*90)
BestOE  = np.array([0]*90)
WorstWT = np.array([0]*90)	
WorstOE = np.array([0]*90)	

Out = open('File.txt','w')

for IdxG in Genes:
	if IdxG[0] == 'Y':
		Genes[IdxG].ProfilingDensity(ProfDensityDict, ProfSamples)

		if Genes[IdxG].CAIBest > 0.9:
			BestPos  = Genes[IdxG].PosBest*3

			if (BestPos - 30 > 0) and (BestPos + 60 < Genes[IdxG].length):
				WT = Genes[IdxG].profdensity[0][BestPos-30:(BestPos + 60)]
				OE = Genes[IdxG].profdensity[2][BestPos-30:(BestPos + 60)]
				if (np.sum(WT) > 0) and (np.sum(OE) > 0):
					BestWT  += WT*1000.0/np.sum(WT)
					BestOE  += OE*1000.0/np.sum(OE)

		if Genes[IdxG].CAIWorst < 0.2:
			WorstPos = Genes[IdxG].PosWorst*3

			if (WorstPos - 30 > 0) and (WorstPos + 60 < Genes[IdxG].length):
				WT = Genes[IdxG].profdensity[0][WorstPos-30:(WorstPos + 60)]
				OE = Genes[IdxG].profdensity[2][WorstPos-30:(WorstPos + 60)]
				if (np.sum(WT) > 0) and (np.sum(OE) > 0):
					WorstWT += WT*1000.0/np.sum(WT)
					WorstOE += OE*1000.0/np.sum(OE)

for Idx in xrange(BestWT.size):
	Out.write(str(BestWT[Idx]) + '\t' + str(BestOE[Idx]) + '\t')
	Out.write(str(WorstWT[Idx]) + '\t' + str(WorstOE[Idx]) + '\n')


