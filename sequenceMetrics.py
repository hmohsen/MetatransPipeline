# By Hussein Mohsen
# Capstone Project in Bioinformatics
# June 2014

import sys

# updates the list of sequence lengths after it parses the FSTA file
def updateList(fileName, list):
	seqsFile = open(fileName)

	line = seqsFile.readline().strip()
	length=0

	while 1:
		if not line:
			list.append(length)
			break

		if line.startswith('>'):
			if length>0:
				list.append(length)
				length=0
		else:
			length=length+len(line)

		line = seqsFile.readline().strip()

	seqsFile.close()

# calculates and print the metrics of the transcripts file (N50, average sequence length, min sequence length, max sequence length, median sequence length, and number of sequence)
def metrics(seqsList):
	total=sum(seqsList)

	#sort in decreasing order
	seqsList.sort(reverse=True)
	
	sumSoFar=0

	print "================================================"
	#Number of sequences
	print "Number of sequences: "+str(len(seqsList))

	# N50 and L50
	# DEFINITIONS ACCORDING TO FAMOUS ASSEMBLATHON PAPER:
	# N50: Sort sequences in decreasing order, loop over sequences and sum them. N50 is the length of the sequence whose addition to the sum leads the summed value to reach one-half of total length of all sequences.
	# L50: Number of sequences in the assembly whose length is greater than or equal to N50
	for i in range(len(seqsList)):
		sumSoFar=sumSoFar+seqsList[i]

		if sumSoFar>=total/2:
			N50=seqsList[i]
			print "N50 value of sequences: "+str(N50)
			
			#all sequences before sequence i are of equal or greater lenth of it since list is sorted in decreasing order
			L50=i+1

			# now check if there are sequences that have equal size to sequence i before it
			for j in range(i+1,len(seqsList)):
				if seqsList[j]==N50:
					L50=L50+1
				elif seqsList[j]<N50:
					break #breaks from inner loop


			print "L50 value of sequences: "+str(L50)
			break #breaks from outer loop

	# Average sequence length
	print "Average sequence length: "+str(total/len(seqsList))

	# Median sequence length (checks below if length of list is even or odd to calculate median accordingly)
	median=seqsList[len(seqsList)/2]

	if len(seqsList)%2==0:
		print "Median sequence length: "+str(median)
	else:
		print "Median sequence length: "+str((median+seqsList[(len(seqsList)/2)-1])/2)

	# Maximum and minimum sequence length
	print "Minimum sequence length: "+str(seqsList[-1])
	print "Maximum sequence length: "+str(seqsList[0])
	print "================================================"
	print "NOTE THAT N50 and L50 ARE CALCULATED ABOVE ACCORDING TO THEIR DEFINITIONS IN FAMOUS ASSEMBLATHON PAPER:\n(http://assemblathon.org/post/44431933387/assemblathon-2-basic-assembly-metrics)"
	print "- N50: Sort contigs in decreasing order, loop over contigs and sum them. N50 is the length of the contig whose addition to the sum leads the summed value to reach one-half of total length of all contigs."
	print "- L50: Number of contigs in the assembly whose length is greater than or equal to N50."



# main method.
def main():
	if __name__ == "__main__":
		seqsList=[]

		for i in range(1, len(sys.argv)):
			seqFileName=str(sys.argv[i])
			updateList(seqFileName, seqsList)

		metrics(seqsList)
		print "Done."

main()
