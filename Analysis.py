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

#ProfilingSamples = ['../WIGs/FC278-L6-P1-ACTTGA-WT','../WIGs/FC278-L5-P1-TGACCA-DhKO','../WIGs/FC278-L5-P1-CAGATC-Old-DhOE']
#RNASeqSamples    = ['../WIGs/HBGN1ADXX-2-ACTTGA-WT','../WIGs/HBGN1ADXX-1-TGACCA-DhKO','../WIGs/HBGN1ADXX-1-CAGATC-DhOE']

#ProfilingDensityDict = {}
#RNASeqDensityDict    = {}   

#for Idx in ProfilingSamples:
	#print Idx
	#ProfilingDensityDict[Idx] = LoadDict(Idx, LengthDict)

#for Idx in RNASeqSamples:
	#print Idx
	#RNASeqDensityDict[Idx]    = LoadDict(Idx, LengthDict)		

# Generate a list of all the gene sequences

File  = open('1-CDS-List/Coding.txt')
Genes = {}
    
for line in File:
	Name = line.split('\t')[1]
	Genes[Name] = Gene(line, GODict, Expression, HalfLives)
	Genes[Name].DefineSeq(Genome[Genes[Name].chr])
	Genes[Name].CalculateCAI(CAIDict)
	#Genes[Name].ProfilingDensity(ProfilingDensityDict, ProfilingSamples)
	#Genes[Name].RNASeqDensity(RNASeqDensityDict, RNASeqSamples)

File.close()

High = ['YMR230W','YBR181C','YKL056C','YJL190C','YML026C','YOR063W','YPL090C','YEL034W','YHL015W','YLR325C','YNL067W','YDL191W','YLL024C','YOL121C','YDR450W','YGL135W','YGL123W','YNL302C','YOR293W','YGR254W','YDL083C','YGL147C','YOL086C','YKL152C','YHR010W','YIL052C','YDR385W','YOR369C','YGL103W','YGR192C','YPR080W','YPL131W','YGL030W','YIL018W','YLR029C','YBR118W','YJL138C','YKL180W']
All  = ['YLR355C','YJR123W','YDR447C','YHR174W','YMR230W','YBR181C','YKL056C','YJL190C','YML024W','YML026C','YLR167W','YOR063W','YPL090C','YOL040C','YLL045C','YBR189W','YDR382W','YEL034W','YHL015W','YLR325C','YLR134W','YHR141C','YLR075W','YKL060C','YOL039W','YPR102C','YCR012W','YOR133W','YNL067W','YGL031C','YDL191W','YAL003W','YBL092W','YLR340W','YLL024C','YOL121C','YML063W','YDR450W','YLR061W','YBR031W','YPR043W','YOL127W','YPR132W','YDR012W','YPL249C-A','YGL076C','YGL135W','YNL209W','YJL189W','YGL123W','YNL302C','YOR293W','YNL069C','YMR142C','YGR254W','YDL083C','YDR418W','YGL189C','YIL069C','YGL147C','YOL086C','YKL152C','YPL198W','YHR010W','YDR064W','YIL052C','YHL033C','YOR096W','YDR385W','YGR148C','YAL038W','YOR369C','YGL103W','YGR192C','YJR009C','YJL052W','YPR080W','YLR249W','YPL131W','YPL220W','YKR059W','YMR116C','YGL030W','YDR050C','YLR110C','YNL178W','YCR031C','YLR044C','YER074W','YIL018W','YGR118W','YOL120C','YLR029C','YDR500C','YBR118W','YGL008C','YLR048W','YJL138C','YDL229W','YKL180W']	
	
for Idx in All:
	if Idx  in High:
		print Idx
			
#FileOut = open('Test.txt','w')

#FileOut.write('Gene\tCAI\tWTProf\tWTRNA\tKOProf\tKORNA\tOEProf\tOERNA\n')

#for IdxG in Genes:
	#FlagCheck = 1.0
	#for Idx in xrange(0,3):
		#FlagCheck = FlagCheck*Genes[IdxG].profRPKM[Idx]*Genes[IdxG].mRNARPKM[Idx]
	#if FlagCheck != 0.0:
		#if IdxG[0] == 'Y' or IdxG[0] == 'Q':
			#if Genes[IdxG].CAI < 0.3:
				#FileOut.write(IdxG + '\t0.3\t')
			#else:
				#FileOut.write(IdxG + '\t' + '{0:.1f}'.format(Genes[IdxG].CAI) + '\t')
			#for Idx in xrange(0,3):
				#FileOut.write(str(Genes[IdxG].profRPKM[Idx]) + '\t' + str(Genes[IdxG].mRNARPKM[Idx]) + '\t')
		
			#FileOut.write('\n')
			
#print time.time() - Start

	
