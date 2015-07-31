import subprocess as sp
import time
import pickle

from itertools import izip

PHREDDict = {
	'!':9.999999e-01, '"':7.943282e-01, '#':6.309573e-01, '$':5.011872e-01, '%':3.981072e-01,
	'&':3.162278e-01, '\'':2.511886e-01, '(':1.995262e-01, ')':1.584893e-01, '*':1.258925e-01,
	'+':1.000000e-01, ',':7.943282e-02, '-':6.309573e-02, '.':5.011872e-02, '/':3.981072e-02,
	'0':3.162278e-02, '1':2.511886e-02, '2':1.995262e-02, '3':1.584893e-02, '4':1.258925e-02,
	'5':1.000000e-02, '6':7.943282e-03, '7':6.309573e-03, '8':5.011872e-03, '9':3.981072e-03,
	':':3.162278e-03, ';':2.511886e-03, '<':1.995262e-03, '=':1.584893e-03, '>':1.258925e-03,
	'?':1.000000e-03, '@':7.943282e-04, 'A':6.309573e-04, 'B':5.011872e-04, 'C':3.981072e-04,
	'D':3.162278e-04, 'E':2.511886e-04, 'F':1.995262e-04, 'G':1.584893e-04, 'H':1.258925e-04,
	'I':1.000000e-04, 'J':7.943282e-05
	}
	
NameDict = {
	  'Scchr01':'chrI',   'Scchr02':'chrII', 'Scchr03':'chrIII',
	 'Scchr04':'chrIV',    'Scchr05':'chrV',  'Scchr06':'chrVI',
	'Scchr07':'chrVII', 'Scchr08':'chrVIII',  'Scchr09':'chrIX',
	  'Scchr10':'chrX',   'Scchr11':'chrXI', 'Scchr12':'chrXII',
   'Scchr13':'chrXIII',  'Scchr14':'chrXIV',  'Scchr15':'chrXV',
	'Scchr16':'chrXVI',     'Scmito':'chrM'
}

LengthDict = {
	  'chrI':230218,   'chrII':813184,  'chrIII':316620, 'chrIV':1531933,   'chrV':576874,  
	 'chrVI':270161, 'chrVII':1090940, 'chrVIII':562643,  'chrIX':439888,   'chrX':745751, 
     'chrXI':666816, 'chrXII':1078177, 'chrXIII':924431, 'chrXIV':784333, 'chrXV':1091291,
    'chrXVI':948066,     'chrM':85779
}

Files    = ['HBGN1ADXX-1-CAGATC-DhOE','HBGN1ADXX-1-TGACCA-DhKO','HBGN1ADXX-2-ACTTGA-WT','HBGN1ADXX-2-TGACCA-T+']

for Idx in Files:
	
	Start = time.time()
	
	# Perform cutadapt Adapter removal
	
	CutAdapt = ["cutadapt","-a","CTGTAGGCACCATCAATAGATCGGAA","-o","2-CutAdapt/","1-Unprocessed/"]
		
	CutAdapt[4] += Idx + '.fastq' 
	CutAdapt[5] += Idx + '.fastq.gz'
	Process = sp.Popen(CutAdapt, stdout=sp.PIPE, stderr=sp.PIPE)
	
	# Take the output of the CutAdapt call and print it to file!
	
	FileOut = open('2-CutAdapt/Reports/' + Idx + '.txt','w')
	FileOut.write(Process.communicate()[0])
	FileOut.close()
	
	print Idx, ' - CutAdapt - ', time.time() - Start
	
	# Get number of lines and divide by four to get the number of individual entries
	
	Length  = int(sp.Popen(["wc","-l","2-CutAdapt/" + Idx + ".fastq"], stdout=sp.PIPE).communicate()[0].split(' ')[-2])/4
	
	File    = open('2-CutAdapt/' + Idx + '.fastq')
	FileOut = open('3-Filtered/' + Idx + '.fastq','w')
	Report  = open('3-Filtered/CutList/' + Idx + '-Stripped.txt','w') 
	
	# Calculate the quality score for each line and trim everything with less than a 99.5% confidence in sequence
	#
	# The confidence in the read is defined as the product of (1 - the quality scores of the individual nucleotides)
	
	for IdxN in range(0,Length):
		Identifier  = File.readline().rstrip('\n')
		Sequence    = File.readline().rstrip('\n')
		QIdentifier = File.readline().rstrip('\n')
		PHRED       = File.readline().rstrip('\n')
		Score       = 1.0
		Length      = len(PHRED)
		
		if Length < 4:
			Report.write('{0:9d}'.format(IdxN) + '\t' + str(Length) + '\t' + Sequence + '\n')
			pass
		else:
			for IdxL in range(0,Length):
				Score   = Score*(1 - PHREDDict[PHRED[IdxL]])
			
			if (Score > 0.995):
				FileOut.write(Identifier + '\n' + Sequence + '\n' + QIdentifier + '\n' + PHRED + '\n')
			else:
				Report.write('{0:9d}'.format(IdxN) + '\t' + str(Length) + '\t' + Sequence + '\n')
			
	File.close()
	FileOut.close()
	Report.close()
	
	print Idx, ' - QualFilt - ', time.time() - Start
	
	sp.call(['gzip','2-CutAdapt/' + Idx + '.fastq'])
	
	# Use Bowtie to align to the non-coding RNAs/tRNAs...
	#
	# The bowtie settings being used are...
	#
	
	Input   = '3-Filtered/' + Idx + '.fastq'
	Output  = '4-NC-Subtracted/SAM/' + Idx + '.sam'
	Bowtie  = ['0-Tools/0-bowtie/bowtie','-Sv','2','--best','0-Tools/2-Indexes/Yeast-Noncoding/Yeast-Noncoding',Input,Output]
	Process = sp.Popen(Bowtie, stdout=sp.PIPE, stderr=sp.PIPE)
	
	FileOut = open('4-NC-Subtracted/Reports/' + Idx + '.txt','w')
	FileOut.write(Process.communicate()[1])
	FileOut.close()
	
	print Idx, ' - BowtieNC - ', time.time() - Start
	
	# ...and then iterate through the FASTQ file to remove the lines that map to nc RNAs/tRNAs!
	
	Length  = int(sp.Popen(["wc","-l","3-Filtered/" + Idx + ".fastq"], stdout=sp.PIPE).communicate()[0].split(' ')[-2])/4
	
	Filtered   = open('3-Filtered/' + Idx + '.fastq')
	Matches    = open('4-NC-Subtracted/SAM/' + Idx + '.sam')
	Subtracted = open('4-NC-Subtracted/' + Idx + '.fastq','w')
	Report     = open('4-NC-Subtracted/CutList/' + Idx + '-Stripped.txt','w')
	
	Count      = 0
	
	for Burn in range(0,415):
		Matches.readline()
		
	for Iter in range(0,Length):
		Count   += 1
		
		Identifier  = Filtered.readline().rstrip('\n')
		Sequence    = Filtered.readline().rstrip('\n')
		QIdentifier = Filtered.readline().rstrip('\n')
		PHRED       = Filtered.readline().rstrip('\n')
		
		SplitF      = Identifier.split(' ')
		SplitM      = Matches.readline().split('\t')
		
		if (SplitM[1] == '4'):
			Subtracted.write(Identifier + '\n' + Sequence + '\n' + QIdentifier + '\n' + PHRED + '\n')
		else:
			Report.write('{0:9d}'.format(Count)+ '\t' + str(len(Sequence)) + '\t' + Sequence + '\n')
			
	Filtered.close()
	Matches.close()
	Subtracted.close()
	Report.close()
	
	print Idx, ' - FilterNC - ', time.time() - Start
	
	# Use Bowtie to align to the genome!
	
	Input   = '4-NC-Subtracted/' + Idx + '.fastq'
	Output  = '5-Aligned/' + Idx + '.sam'
	Bowtie  = ['0-Tools/0-bowtie/bowtie','-Syam','1','--best','--strata','0-Tools/2-Indexes/SCerevisiae/s_cerevisiae', Input, Output]
	Process = sp.Popen(Bowtie, stdout=sp.PIPE, stderr=sp.PIPE)
	
	FileOut = open('5-Aligned/Reports/' + Idx + '.txt','w')
	FileOut.write(Process.communicate()[1])
	FileOut.close()
	
	print Idx, ' - BowtieGe - ', time.time() - Start
	
	sp.call(['gzip','3-Filtered/' + Idx + '.fastq'])
	sp.call(['gzip','4-NC-Subtracted/' + Idx + '.fastq'])
	sp.call(['gzip','4-NC-Subtracted/SAM/' + Idx + '.sam'])

	File       = open('5-Aligned/' + Idx + '.sam')
	
	PlusDict   = {}
	MinusDict  = {}
	OffsetDict = {
		21:11, 22:11, 23:11, 24:12, 25:12, 26:12,
		27:13, 28:13, 29:13, 30:14, 31:14, 32:14,
		33:15, 34:15, 36:15
	}
	
	Count = 0.0
	
	for IdxC in LengthDict:
		PlusDict[IdxC]  = [0]*LengthDict[IdxC]
		MinusDict[IdxC] = [0]*LengthDict[IdxC]

	for line in File:
		if line[0] != '@':
			Split = line.split('\t')
			
			if len(Split[9]) in OffsetDict:
				Offset = OffsetDict[len(Split[9])] - 1
			else:
				Offset = -1
			
			if Split[1] == '16':
				MinusDict[NameDict[Split[2]]][int(Split[3]) + Offset] += 1
				Count += 1
			elif Split[1] == '0':
				PlusDict[NameDict[Split[2]]][int(Split[3]) + Offset] += 1
				Count += 1
			else:
				pass
			
	File.close()
	
	FileOutP   = open('6-WIGs/' + Idx + '-P.wig','w')
	FileOutM   = open('6-WIGs/' + Idx + '-M.wig','w')
	
	for IdxC in PlusDict:
		FileOutP.write('track name=tracklabel viewLimits=-5:5 color=79,159,36\n')
		FileOutM.write('track name=tracklabel viewLimits=-5:5 color=79,159,36\n')
		
		FileOutP.write('fixedStep  chrom=' + IdxC + '  start=1  step=1\n')
		FileOutM.write('fixedStep  chrom=' + IdxC + '  start=1  step=1\n')
		
		for IdxN in range(0, len(PlusDict[IdxC])):
			FileOutP.write(str(PlusDict[IdxC][IdxN]/Count*1000000) + '\n')
			FileOutM.write(str(MinusDict[IdxC][IdxN]/Count*1000000) + '\n')
	
	FileOutP.close()	
	FileOutM.close()

	print Idx, ' - WIGsMade - ', time.time() - Start