from Gene import *
from SequenceTools import *

import pickle

# Load the reference genome

PickleNames   = open('1-CDS-List/NameConvert.pckl')
PickleGenome  = open('2-Ref-Genome/Genome.pckl')
PickleGODict  = open('4-Ontologies/GOSlim.pckl')
PickleExpress = open('5-Expression/Expression.pckl')

Names         = pickle.load(PickleNames)
Genome        = pickle.load(PickleGenome)
GODict        = pickle.load(PickleGODict)
Expression    = pickle.load(PickleExpress)
 
PickleNames.close()
PickleGenome.close()
PickleGODict.close()
PickleExpress.close()

# Generate a list of all the gene sequences

File  = open('1-CDS-List/Coding.txt')
Genes = {}
    
for line in File:
	Name = line.split('\t')[1]
	Genes[Name] = Gene(line)
	Genes[Name].DefineSeq(Genome[Genes[Name].chr])
	
	# Ontology assign
	
	if Name in GODict:
		Genes[Name].DefineOntologies(GODict[Name])
	else:
		Genes[Name].DefineOntologies([])
		
	# Expression assign
	
	if Name in Expression:
		Genes[Name].DefineExpression(Expression[Name])
	else:
		Genes[Name].DefineExpression(0.0)

File.close()
    
# Calculate frequency tables!

List = []

for IdxG in Genes:
	if Genes[IdxG].expression > 100000:
		List.append(IdxG)

CAIGenes = FrequencyTable(Genes,List)
AllGenes = FrequencyTable(Genes)

for Idx in sorted(CAIGenes.keys()):
	print Idx, CAIGenes[Idx][0], '{:6.3f}'.format(CAIGenes[Idx][1]/AllGenes[Idx][1])

print CAIGenes['AGA'][1]
print AllGenes['AGA'][1]		