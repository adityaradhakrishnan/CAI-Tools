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