from Gene import Gene

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
	
AbundanceBuild()