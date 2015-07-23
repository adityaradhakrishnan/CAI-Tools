from Gene import Gene
import SequenceTools
import pickle

# Load the reference genome

PickleNames   = open('1-CDS-List/NameConvert.pckl')
PickleGenome  = open('2-Ref-Genome/Genome.pckl')
PickleGODict  = open('4-Ontologies/GOSlim.pckl')
PickleExpress = open('5-Expression/Expression.pckl')
PickleHL      = open('6-RNAHalfLives/HalfLives.pckl')

Names         = pickle.load(PickleNames)
Genome        = pickle.load(PickleGenome)
GODict        = pickle.load(PickleGODict)
Expression    = pickle.load(PickleExpress)
HalfLives     = pickle.load(PickleHL)
 
PickleNames.close()
PickleGenome.close()
PickleGODict.close()
PickleExpress.close()
PickleHL.close()

# Generate a list of all the gene sequences

File  = open('1-CDS-List/Coding.txt')
Genes = {}
    
for line in File:
	Name = line.split('\t')[1]
	Genes[Name] = Gene(line, GODict, Expression, HalfLives)
	Genes[Name].DefineSeq(Genome[Genes[Name].chr])
	Genes[Name].UWGenome()
	
File.close()

File = open('Analysis/UWBoundary-RNA.tsv','w')
File.write('Gene\tPolyAHL\tHL\tUWCounts\tUWFrac\n')

for IdxG in Genes:
	if IdxG in HalfLives:
		File.write(IdxG + '\t' + str(Genes[IdxG].halflife[0]) + '\t' + str(Genes[IdxG].halflife[1]) + '\t' + str(Genes[IdxG].uwcount) + '\t' + str(Genes[IdxG].uwfrac) + '\n')
			
File.close()

File = open('Analysis/UWBoundary-Protein.tsv','w')
File.write('Gene\tExpression\tUWCounts\tUWFrac\n')

for IdxG in Genes:
	if IdxG in Expression:
		File.write(IdxG + '\t' + str(Genes[IdxG].expression) + '\t' + str(Genes[IdxG].uwcount) + '\t' + str(Genes[IdxG].uwfrac) + '\n')
			
File.close()