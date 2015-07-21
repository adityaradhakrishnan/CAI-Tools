def FrequencyTable(*arg):
    
	AACodon = {
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
		'TAC':'Y', 'TAT':'Y', 'TAA':'-', 'TAG':'-',
		'TGC':'C', 'TGT':'C', 'TGA':'-', 'TGG':'W',
		}
    
	CodonTable = {}
	AACounts   = {}
	Count      = 0.0
	
	if len(arg) > 1:
		GeneSet = arg[1]
	else:
		GeneSet = list(arg[0].keys())

	for IdxG in GeneSet:
		for IdxC in arg[0][IdxG].codonlist:
			Count += 1
			
			if AACodon[IdxC] in AACounts:
				AACounts[AACodon[IdxC]] += 1
            		else:
				AACounts[AACodon[IdxC]] = 1.0
			if IdxC in CodonTable:
				CodonTable[IdxC] += 1
			else:
				CodonTable[IdxC]  = 1
	
	Frequencies = {}

	for IdxC in CodonTable:
		Frequencies[IdxC] = [AACodon[IdxC], CodonTable[IdxC]/AACounts[AACodon[IdxC]], CodonTable[IdxC]/Count*1000, CodonTable[IdxC]]

	return Frequencies
	
