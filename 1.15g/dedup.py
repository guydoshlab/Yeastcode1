from Bio import SeqIO
def comparefastqs(fastq1,fastq2):

	f_1=open(fastq1)
	reads1=SeqIO.parse(f_1,"fastq")
	f_2=open(fastq2)
	reads2=SeqIO.parse(f_2,"fastq")
	
	dict1={}
	dict2={}
	tot=0
	tot1=0
	tot2=0
	for read in reads1:
		dict1[str(read.seq)]=1
		tot1+=1
	print "Starting comparison"
	for read in reads2:
		tot2+=1
		
		#if dict2.has_key(str(read.seq)):
		#	print "duplicate read found in dict2: "+str(read.seq)
			
		dict2[str(read.seq)]=1
		
		if dict1.has_key(str(read.seq)):
			tot+=1
			continue
		else:
			print "Doesn't have read: "+str(read.seq)		
			print tot
	print "Total reads in input1: "+str(tot1)
	print "Total reads in input2: "+str(tot2)
	print "Total reads in in 2 not in 1: "+str(tot2-tot)
	f_1.close()
	f_2.close()
	return
	
def removedups(reads1):
	dict1={}
	totdup=0
	tot1=0
	for read in reads1:
		if dict1.has_key(str(read.seq)):
			totdup+=1
			continue
		else:
			dict1[str(read.seq)]=1
			tot1+=1
			yield read
	print "Tot duplicate reads: "+str(totdup)
	print "Tot reads retained: "+str(tot1)
	
		
def writeuniques(fastq_in,fastq_out):
	readsin=SeqIO.parse(fastq_in,"fastq")
	de_jp_reads=removedups(readsin)
	try:
		SeqIO.write(de_jp_reads,fastq_out,"fastq")
	except TypeError:
		print "No duplicates found"
	
	