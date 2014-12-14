# By Hussein Mohsen
# Capstone Project in Bioinformatics
# June-December 2014

import sys
import os

# updates the list of sequence lengths after it parses the FASTA file
def blastGOLDSeqs(fileName, dbDirectory, resultsDirectory, id):
	
	path=""
	
	# create a directory for BLAST results
	path= resultsDirectory+"blastAllResults"+str(id)+"/"

	if not os.path.exists(path):
    		os.makedirs(path)

	seqsFile = open(fileName)

	line = seqsFile.readline().strip()

	seq=''
	counter=1

	while 1:
		if not line:
			exampleFile = open("example.fa", "w")
			exampleFile.write(">seq "+str(counter)+"\n"+seq)
			exampleFile.close()
			
			os.system("blastn -query example.fa -db "+dbDirectory+" -evalue 1e-30 -out "+path+"results"+str(counter)+".out")
			break

		if line.startswith('>'):
			if seq!='':
				exampleFile = open("example.fa", "w")
				exampleFile.write(">seq "+str(counter)+"\n"+seq)
				exampleFile.close()
				
				os.system("blastn -query example.fa -db "+dbDirectory+" -evalue 1e-10 -out "+path+"results"+str(counter)+".out")
				counter=counter+1
				seq=''
		else:
			seq=seq+line

		line = seqsFile.readline().strip()

	seqsFile.close()

# parses a .out file that contains BLAST results
# these files are generated bu blastAll() above
# returns the name and e-value of top hit
def parseOutFile(fileName):

	noHits=0

	outFile = open(fileName, "r")

	line = outFile.readline().strip()

	# initialize lists
	names=[]
	EValues=[]

	while 1:
		if "No hits found" in line:
			noHits=noHits+1 
			break

		if "Matrix:" in line:
			break

		if '>' in line:
			names.append(line[1:].strip())

			line = outFile.readline().strip()
			line = outFile.readline().strip()

			# this will contain the e-value (revise BLAST output files)
			line = outFile.readline().strip()

			EValues.append(float(line[line.index("Expect = ")+8:]))

		line = outFile.readline().strip()

	outFile.close()

	return names, EValues, noHits
	

# generates statistics about BLAST results obtained by blastAll() above
def stats(directory, threshold, transDict):

	total=0
	totalNoHits=0

	files=os.listdir(directory)
	
	for file in files:
		if file.endswith('.out'):
			# total number of transcripts
			total=total+1

			# gets the transcript on first hit in BLAST results
			transNames, eVals, noHits =parseOutFile(directory+file)
			
			totalNoHits=totalNoHits+noHits

			for i in range(0, len(transNames)):
				transName=transNames[i]
				eVal=eVals[i]

				if transName != 'NO_HITS' and eVal<threshold:
					if transName  in transDict:
	  					transDict[transName]=transDict[transName]+1
					else:
						transDict[transName]=1

	return total, totalNoHits

# gets transcripts and the nodes that constitute each
# returns them in a dictionary: key is transcript full name and value is a list of nodes
def getTranscriptsNodes(fileName):

	allTransAndTheirNodes={}

	seqsFile = open(fileName)

	line = seqsFile.readline().strip()

	# parse Nodes to reach transcripts
	while 1:
		if not line:
			# checking for double not lines to parse the Oases contig-order file perfectly
			line = seqsFile.readline().strip()
			
			if not line:
				break

		if line.startswith('>') and 'Transcript' in line:
			# add this transcript name a sa new key and initialize the node list that will serve as its value
			transcName=line[1:]
			allTransAndTheirNodes[transcName]=[]

			# get the line that contains the description of the transcript
			line = seqsFile.readline().strip()

			allTransAndTheirNodes[transcName].append(int(line[0:line.find(':')]))
			
			index = line.find('->')
			
			while index != -1:
				line=line[index+2:]
				allTransAndTheirNodes[transcName].append(int(line[0:line.find(':')]))
				index = line.find('->')
			
		line = seqsFile.readline().strip()

	seqsFile.close()

	# print allTransAndTheirNodes

	return allTransAndTheirNodes

# saves graph edges and their multiplicity values in a 2D dictionary called multipl
def getGraph(graphFileName):

	# dictionary of multiplicity of graph arcs and another for arcs between nodes (list representation of graph)
	multipl={}

	# flag that helps in parsing the file	
	detectedArcs=False

	graphFile = open(graphFileName)

	line = graphFile.readline().strip()

	while 1:
		if not line or (detectedArcs==True and 'ARC' not in line):
			break

		if 'ARC' in line:
			list = line.split('\t')

			fromNode= abs(int(list[1]))
			toNode= abs(int(list[2]))
			multiplicity= int(list[3])

			multipl[fromNode, toNode]=multiplicity
			multipl[toNode, fromNode]=multiplicity

			detectedArcs=True

		line = graphFile.readline()

	#print multipl

	return multipl

# exports to a new file a pair that corresponds to each transcript that got BLAST hits
# this pair includes total multiplicity and blast hits
# parameter maxBlastHits can be know from global statistics file
# parameter exportNodesPerTrans is boolean: if True then it exports TransNodesHisto.txt in output directory (just simply the distribution of nodes per transcript in transcripts that got BLAST hits)
def exportLociHitsAndScores(outputDirectory, transDict, allTransAndTheirNodes, multipl, maxBlastHits, exportNodesPerTrans):

	total=0.0

	# number of transcripts with 1 node only or that got total multiplicity obtained from graph edges
	counter=0

	# used to create the histograms
	# first row contains average multiplicity values for transcripts that got 1 to N/5 hits, second for those with (N/5)+1 to 2N/5 hits and so on
	histos=[[],[],[],[],[]]

	for trans in transDict:
		try:
			temp=allTransAndTheirNodes[trans]

			if len(temp)==1 or total>0:
				counter += 1
				
				if total>0:
				
					# average multiplicity for a transcript
					average = total/len(allTransAndTheirNodes[trans])

					#print str(trans)+" "+str(total)+" "+str(transDict[trans])
					if (transDict[trans]-1)/5 >= 5:
						histos[4].append(average)
					else:
						histos[(transDict[trans]-1)/5].append(average)
			
			# total multiplicity of edges in a transcript that will be used to 
			# calculate average multiplicity based on number of nodes in each transcript
			total=0.0

			for i in range(0, len(temp)-1):
				for j in range(i+1, len(temp)):
					try:
						toNode=temp[i]
						fromNode=temp[j]

						total += multipl[toNode,fromNode]
					except KeyError:
						continue
		
		except KeyError:
			continue
			
	if outputDirectory.endswith('/')==False:
		outputDirectory+='/'

	# create file of histograms data
	# each row corresponds to a histogram
	HistosFile=open(outputDirectory+"allHistos.txt", "w")

	print "Working on and writing in statistics file allHistos.txt..."

	for i in range(0, len(histos)):
		for j in range(0, len(histos[i])):
			if j==len(histos[i])-1:
				HistosFile.write(str(histos[i][j])+'\n')
			else:
				HistosFile.write(str(histos[i][j])+',')

	HistosFile.close()

	#print histos
	#print counter

	if exportNodesPerTrans==True:

		print "Working on and writing in statistics file TransNodesHisto.txt..."

		tempString=''
		
		for val in allTransAndTheirNodes.keys():
			try:
				# we care about transcripts that got BLAST hits
				if transDict[val]>=1:
					tempString+= str(len(allTransAndTheirNodes[val]))+', '
			except KeyError:
				continue

		tempString= tempString[0: len(tempString)-2]

		# create file of histograms data for the distribution of number of nodes per transcript in transcripts that got BLAST hits
		NodesPerTransHistosFile=open(outputDirectory+"TransNodesHisto.txt", "w")

		NodesPerTransHistosFile.write(tempString+'\n')

		NodesPerTransHistosFile.close()

# It exports parameter-outputDirectory/LociTransHisto.txt in output directory (just simply the distribution of number loci sizes (number of transcripts per locus) in parameter-directory/contig-ordering.txt
def exportTransPerLocus(contigOrderingFileName, outputDirectory):

	print "Working on and writing in statistics file LociTransHisto.txt..."

	# created to read from contig-ordering.txt file
	contigOrderingFile=open(contigOrderingFileName, 'r')

	# create file of histograms data for the number-of-transcripts-per-locus-distribution
	transPerLocusHistosFile=open(outputDirectory+"LociTransHisto.txt", "w")

	lociCounter=0
	
	# how many loci that generated more than one transcript
	multiTranLoci=0

	# distribution
	transPerLocus = ''
	
	line = contigOrderingFile.readline().strip()

	# parse Nodes to reach transcripts
	while 1:
		if not line:
			# checking for double not lines to parse the Oases contig-order file perfectly
			line = contigOrderingFile.readline().strip()

			if not line:
				break

		if line.startswith('>') and '_Transcript' in line:
			
			locusNumber = int(line[line.index('_')+1: line.index('_Transcript')])

			# detected transcripts from new locus
			if locusNumber > lociCounter:
				lociCounter = locusNumber
			
				index=line.index('/')

				# how many transcripts generated from this locus
				num = int(line[index+1:line.index('_Confidence')])

				if num>1:
					multiTranLoci+=1
			
				transPerLocus += str(num)+', '
		
		line = contigOrderingFile.readline().strip()

	contigOrderingFile.close()

	transPerLocus=transPerLocus[0: len(transPerLocus)-2]

	transPerLocusHistosFile.write(transPerLocus+'\n')

	print str(multiTranLoci)+' loci with more than 1 transcript.'

# It exports a file parameter-outputDirectory/LociTransBlastHits.txt that describes the # of blast hits that each locus transcripts got
def exportLociTransHits(contigOrderingFileName, outputDirectory, transDict):

	print "Working on and writing in statistics file LociTransBlastHits.txt..."

	# created to read from contig-ordering.txt file
	contigOrderingFile=open(contigOrderingFileName, 'r')

	# create file of  for the number-of-blast-hits-per-transcript-in-each-locus-distribution
	transPerLocusBlastHitsFile=open(outputDirectory+"LociTransBlastHits.txt", "w")

	lociCounter=0

	# distribution per each locus
	hitsPerLocusTrans = ''
	
	line = contigOrderingFile.readline().strip()

	# parse Nodes to reach transcripts
	while 1:
		if not line:
			# checking for double not lines to parse the Oases contig-order file perfectly
			line = contigOrderingFile.readline().strip()

			if not line:
				break

		if line.startswith('>') and '_Transcript' in line:
			
			locusNumber = int(line[line.index('_')+1: line.index('_Transcript')])

			# detected transcripts from new locus
			if locusNumber > lociCounter:
				if lociCounter!=0:
					transPerLocusBlastHitsFile.write('Locus '+str(lociCounter)+' '+hitsPerLocusTrans+'\n')
					hitsPerLocusTrans =''

				lociCounter=locusNumber
			
			transName=line[1:]
			#print transName

			if transName in transDict.keys():
				hitsPerLocusTrans+=(str(transDict[transName])+' ')
			else:
				hitsPerLocusTrans+='0 '

		line = contigOrderingFile.readline().strip()

	transPerLocusBlastHitsFile.write('Locus '+str(locusNumber)+' '+hitsPerLocusTrans+'\n')
	
	contigOrderingFile.close()

# It exports files parameter-outputDirectory/1TransZeroHitsLoci.txt and parameter-outputDirectory/1TransPositiveHitsLoci.txt
# which are the #s of Loci with 1 Oases transcripts that got 0 or >=1 BLAST hits, respectively
# This can be used by method export1TransLociScores() in EM.py
def export1TransLoci(contigOrderingFileName, outputDirectory, transDict):

	print "Working on 1TransLoci files..."

	# created to read from contig-ordering.txt file
	contigOrderingFile=open(contigOrderingFileName, 'r')

	# create files described
	zeroFile=open(outputDirectory+"1TransZeroHitsLoci.txt", "w")
	positFile=open(outputDirectory+"1TransPositiveHitsLoci.txt", "w")

	lociCounter=0

	# distribution per each locus
	hitsPerLocusTrans = ''
	
	line = contigOrderingFile.readline().strip()

	# parse Nodes to reach transcripts
	while 1:
		if not line:
			# checking for double not lines to parse the Oases contig-order file perfectly
			line = contigOrderingFile.readline().strip()

			if not line:
				break

		if line.startswith('>') and '_Transcript_1/1_' in line:
			
			locusNumber = int(line[line.index('_')+1: line.index('_Transcript')])

			# detected transcripts from new locus
			if locusNumber > lociCounter:
				if lociCounter!=0:
					if ' 0' in hitsPerLocusTrans:
						zeroFile.write(str(lociCounter)+' ')
					else:
						positFile.write(str(lociCounter)+' ')
					
					hitsPerLocusTrans =''

				lociCounter=locusNumber
			
			transName=line[1:]
			#print transName

			if transName in transDict.keys():
				hitsPerLocusTrans+=(str(transDict[transName])+' ')
			else:
				hitsPerLocusTrans+=' 0'

		line = contigOrderingFile.readline().strip()

	contigOrderingFile.close()
	zeroFile.close()
	positFile.close()

# main method.
def main():
	if __name__ == "__main__":

		funcsType= 1

		# total number of GOLD transcripts		
		total=0

		# total number of GOLD transcripts	that did not get any BLAST results (returned "No Hits found." by BLAST)
		totalNoHits=0

		# dictionary of transcripts
		transDict={}

		start = 1

		# order to determine which functions are needed from this script. Check READ ME file for more details.
		index=sys.argv.index("-f")

		funcsType= int(sys.argv[index+1])

		# order to determine the results directory to fild the transcripts file to be used as BLAST database in blastGoldSeqs()
		index=sys.argv.index("-d")

		resultsDirectory = sys.argv[index+1]

		start=5

		# Statistics (and blastAll if funcsType >2)
		if funcsType >0:
			for i in range(start, len(sys.argv)):
				seqsFileName=str(sys.argv[i])
				
				# blastAll
				if funcsType >2:
					print "Working on BLAST of "+seqsFileName+"..."
					blastGOLDSeqs(seqsFileName, resultsDirectory+"res_mergedAssembly/transcripts.fa", resultsDirectory, (i-start+1))
				
				print "Working on Statistics of "+seqsFileName+"..."
				tempTotal, tempNoHits = stats(resultsDirectory+"blastAllResults"+str(i-start+1)+"/", 1e-30, transDict)
				total=total+tempTotal
				totalNoHits=totalNoHits+tempNoHits

			print "Writing in statistics file blastAllGLOBALStats.txt..."
			# GLOBAL statistics to be saved in running directory #

			# check how many transcripts generated by Oases were matched as first hit exactly once while BLASTing GOLD transcripts (value of singleGuess)
			# check how many transcripts generated by Oases were matched as first hit more than once while BLASTing GOLD transcripts (value of multipleGuess)
			# calculate total number of sequences that got BLAST hits (value of totalSeqs)

			totalHits=0
			singleGuess=0
			multipleGuess=0

			maxHits=0

			for key in transDict.keys():
				totalHits=totalHits+transDict[key]
	
				if(transDict[key]==1):
					singleGuess = singleGuess+1
				elif(transDict[key]>1):
					multipleGuess = multipleGuess+1

					if transDict[key]>maxHits:
						maxHits=transDict[key]

			# Final statistics file that sums up numbers from all local blast results (if there are multiple input files given in command) 
			# blastAllGLOBALStats.txt is created in running directory
			finalStatsFile=open(resultsDirectory+"blastAllGLOBALStats.txt", "w")
			percent = (float(total-totalNoHits)/total)*100
			finalStatsFile.write("Total Accuracy of Oases Assembly: "+str(percent)+"%\n")
			finalStatsFile.write("Total Predicted Errors in Oases Assembly: "+str(100-percent)+"%\n")
			finalStatsFile.write("Total number of Oases sequences that were BLAST-hit more than once: "+str(multipleGuess)+"\n")
			finalStatsFile.write("Maximum number of BLAST-hit one Oases sequence got: "+str(maxHits)+"\n")
	
			finalStatsFile.close() 

			# exports the file that contains how many hit each transcript in each locus got
			exportLociTransHits(resultsDirectory+'res_mergedAssembly/contig-ordering.txt', resultsDirectory+'res_mergedAssembly/', transDict)
			export1TransLoci(resultsDirectory+'res_mergedAssembly/contig-ordering.txt', resultsDirectory+'res_mergedAssembly/', transDict)
		
		# Getting graph and multiplicities to generate data for the histogram
		# that will help us determine significance of low-multiplicity transcripts
		if funcsType >1:
			allTransAndTheirNodes=getTranscriptsNodes(resultsDirectory+'res_mergedAssembly/contig-ordering.txt')
			multipl= getGraph(resultsDirectory+'res_mergedAssembly/LastGraph')
			
			exportLociHitsAndScores(resultsDirectory+'res_mergedAssembly/', transDict, allTransAndTheirNodes, multipl, 26, True)
			exportTransPerLocus(resultsDirectory+'res_mergedAssembly/contig-ordering.txt', resultsDirectory+'res_mergedAssembly/')

			#print allTransAndTheirNodes
			#print multipl

		print "Done."

main()
