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


String = '../Sequencing/2014-12-Dhh1RNASeq/6-WIGs/'
List   = ['HBGN1ADXX-2-ACTTGA-WT','HBGN1ADXX-1-TGACCA-DhKO','HBGN1ADXX-1-CAGATC-DhOE','HBGN1ADXX-2-TGACCA-T+']
#ProfilingSamples = ['../WIGs/Profiling/FC278-L6-P1-ACTTGA-WT','../WIGs/Profiling/FC278-L5-P1-TGACCA-DhKO','../WIGs/Profiling/FC278-L5-P1-CAGATC-Old-DhOE']
RNASeqSamples    = [String + Idx for Idx in List]

#ProfilingDensityDict = {}
RNASeqDensityDict    = {}   

#for Idx in ProfilingSamples:
	#print Idx
	#ProfilingDensityDict[Idx] = LoadDict(Idx, LengthDict)

for Idx in RNASeqSamples:
	print Idx
	RNASeqDensityDict[Idx]    = LoadDict(Idx, LengthDict)		

# Generate a list of all the gene sequences

File  = open('1-CDS-List/Coding.txt')
Genes = {}
    
for line in File:
	Name = line.split('\t')[1]
	Genes[Name] = Gene(line, GODict, Expression, HalfLives)
	Genes[Name].DefineSeq(Genome[Genes[Name].chr])
	Genes[Name].CalculateCAI(CAIDict)
	#Genes[Name].ProfilingDensity(ProfilingDensityDict, ProfilingSamples)
	Genes[Name].RNASeqDensity(RNASeqDensityDict, RNASeqSamples)

FileOut = open('mRNARPKM.txt','w')
		
for IdxG in Genes:
	FileOut.write(IdxG + '\t' + str(Genes[IdxG].length))
	for Idx in Genes[IdxG].mRNARPKM:
		FileOut.write('\t' + str(Idx))

	FileOut.write('\n')
