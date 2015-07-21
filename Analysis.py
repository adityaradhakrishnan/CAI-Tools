from Gene import *
from SequenceTools import *

import pickle

# Load the reference genome

PickleGenome = open('2-Ref-Genome/Genome.pckl')
Genome       = pickle.load(PickleGenome)
PickleGenome.close()

# Generate a list of all the gene sequences

File  = open('1-CDS-List/Coding.txt')

Genes  = {}
    
for line in File:
	Name = line.split('\t')[1]
	Genes[Name] = Gene(line)
	Genes[Name].DefineSeq(Genome[Genes[Name].chr])
    
# Calculate frequency tables!

CodonTable = FrequencyTable(Genes)

for IdxC in sorted(list(CodonTable.keys())):
    print IdxC, CodonTable[IdxC][0], CodonTable[IdxC][1], CodonTable[IdxC][2], CodonTable[IdxC][3]