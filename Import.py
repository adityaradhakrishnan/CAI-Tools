from Gene import Gene

import numpy as np
import pickle

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
	
def NameConversionTable(InFile):
	File     = open(InFile)
	NameDict = {}
	
	for line in File:
		Split = line.rstrip('\n').split(' ')
		
		if Split[4] == '\"\"':
			pass
		else:
			NameDict[Split[4]] = Split[2]
			
	File.close()
	
	File = open(InFile.split('.')[0] + '.pckl','w')
	pickle.dump(NameDict, File)
	File.close()
	
def OntologyBuild():
	GODict = {}
	
	for IdxS in ['Component','Function','Process']:
		File = open('4-Ontologies/' + IdxS + '-Slim.csv')
		for line in File:
			Split = line.rstrip('\n').split('\t')
			for Idx in Split[4].split(','):
				if Idx in GODict:
					GODict[Idx].append(int(Split[0].zfill(7)))
				else:
					GODict[Idx] = [int(Split[0].zfill(7))]
					
		File.close()
	
	File = open('4-Ontologies/GOSlim.pckl','w')
	pickle.dump(GODict,File)
	File.close()

def AbundanceBuild():
	AbundanceDict = {}	
	
	File = open('5-Expression/Expression.txt')
	for line in File:
		Split = line.rstrip('\n').split(' ')
		if Split[1] != '\"\"':
			AbundanceDict[Split[0]] = float(Split[1])
	
	File.close()
			
	File = open('5-Expression/Expression.pckl','w')
	pickle.dump(AbundanceDict,File)
	File.close()
	
def HalfLifeBuild():
	HalfLifeDict = {}
	
	File = open('6-RNAHalfLives/HalfLives.tsv')
	File.readline()
	
	for line in File:
		Split = line.rstrip('\n').split('\t')
		HalfLifeDict[Split[0]] = [float(Split[i]) for i in xrange(2,6)]
	
	File.close()
		
	File = open('6-RNAHalfLives/HalfLives.pckl','w')
	pickle.dump(HalfLifeDict,File)
	File.close()
	
def ClipSeqProcess():
	CAIDict    = pickle.load(open('7-CAI/Gene-CAI.pckl'))
	
	Conditions = ['Dhh1-V1','Dhh1-V2','Sbp1-V1','Sbp1-V2','Lsm1-V1','Lsm1-V2','Pat1-V1','Pat1-V2','WT']
	Gene       = {}
	
	for IdxC in Conditions:
		File = open('9-CLIP-Seq/' + IdxC + '.txt')
		File.readline()
		
		for line in File:
			Split = line.rstrip('\n').split('\t')
			if Split[0] in Gene:
				Gene[Split[0]][IdxC]  = np.array([float(Idx) for Idx in Split[int(Split[1])+6:-1*int(Split[5])]])
			else:
				Gene[Split[0]] = {IdxC: np.array([float(Idx) for Idx in Split[int(Split[1])+6:-1*int(Split[5])]])}
	
		File.close()
	
	FileOut = open('9-CLIP-Seq/RPKM.txt','w')
	
	Counts = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
	
	for Idx in range(0,9):
		IdxC = Conditions[Idx]
		for IdxG in Gene:
			Counts[Idx] += np.sum(Gene[IdxG][IdxC])

	print Counts
				
	for IdxG in Gene:
		Check = 1.0
		for IdxC in Gene[IdxG]:
			Check = Check*np.sum(Gene[IdxG][IdxC])
			
		if Check == 0.0:
			pass
		else:
			if IdxG in CAIDict:
				Out = '0.3'	if (round(CAIDict[IdxG]*10) < 4) else str(CAIDict[IdxG])[0:3]
				FileOut.write(IdxG + '\t' + Out + '\t')
				for Idx in range(0,9):
					IdxC = Conditions[Idx]
					RPM  = np.sum(Gene[IdxG][IdxC])/Counts[Idx]*1000000
					FileOut.write(str(RPM/Gene[IdxG][IdxC].size*1000) + '\t')
			
				FileOut.write('\n')
		
		
	
ClipSeqProcess()