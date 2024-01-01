import csv
from collections import defaultdict
import thread
import operator
import re
import os
import sys
import glob
import os.path
countNuc = 0
list_of_files1 = glob.glob('fasta/fasta_filtered/*fasta')
#list_of_files1 = glob.glob('negati/s3No1Allpiranha9_1_2_21AnnotNR.fasta')
for file_name1 in list_of_files1:
        r1 = open(file_name1,'r')
	file_name1 = file_name1.split("/")
	file_name1 = file_name1[2].split(".fasta")
	file_name1 = file_name1[0]+"motifs4.tsv"	
	outfile = open(os.path.join('fasta/fasta_filtered/motifs/', os.path.basename(file_name1)), 'w')
        print(file_name1)
	line1 = r1.readline().strip("\n")
	#line12 = r1.readline().strip()
	#line13 = r1.readline()
	#line14 = r1.readline()

	countNuc = 0	
	motifs={}
	k=4

	alt_map = {'ins':'0'}
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

	def reverse_complement(kmer):
        	for k,v in alt_map.iteritems():
                	kmer = kmer.replace(k,v)
                	bases = list(kmer)
                	bases = reversed([complement.get(base,base) for base in bases])
                	bases = ''.join(bases)
        	for k,v in alt_map.iteritems():
                	bases = bases.replace(v,k)
                	return bases

	while len(line1) > 0:
		if line1[0] != ">":
			for start in range(len(line1) - k + 1):
				kmer = line1[start:start + k]
				current_count = motifs.get(kmer, 0)
				motifs[kmer] = current_count + 1 
			for s in line1:
				countNuc += 1
        	line1 = r1.readline().strip("\n")
        	#line12 = r1.readline().strip()
        	#line13 = r1.readline()
        	#line14 = r1.readline()
	data=[]


	for kmer in motifs:
		if 'N' not in kmer:
			data.append([kmer, motifs[kmer]])
	data.sort(key=operator.itemgetter(1), reverse=True) 
        #print(countNuc)
	for kmer, ct in data:
                #print(countNuc)
		normalizedCount = str(ct/float(countNuc))
		outfile.write(kmer + '\t' + str(ct) + '\t' + reverse_complement(kmer) + "\t" +normalizedCount +'\n' )


	outfile.close()

"""alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def reverse_complement(kmer):
        for k,v in alt_map.iteritems():
                seq = seq.replace(k,v)
                bases = list(kmer)
                bases = reversed([complement.get(base,base) for base in bases])
                bases = ''.join(bases)
        for k,v in alt_map.iteritems():
                bases = bases.replace(v,k)
                return bases"""


""" barcode = line1.strip().split(" ")[1].split(":")[3][:6]
ifive = line1.strip().split(" ")[1].split(":")[3][-8:]
            if i in line14:
                if ifive == 'AGTTGTCG' or ifive == 'AGATCTCG':
                    if barcode not in outFile:
                        outFile[barcode]=gzip.open(os.path.join(outDir,barcode+tail1),'w')
                        outFile2[barcode]=gzip.open(os.path.join(outDir,barcode+tail2),'w')
           `        line12 = r1.readline()
                    line13 = r1.readline()
                    line14 = r1.readline()
                    line22 = r2.readline()
                    line23 = r2.readline()
                    line24 = r2.readline()
                    if len(line12) == len(line14) and len(line22)==len(line24):
                        outFile[barcode].write(line1)
                        outFile[barcode].write(line12)
                        outFile[barcode].write(line13)
                        outFile[barcode].write(line14)
                        outFile2[barcode].write(line2)
                        outFile2[barcode].write(line22)
                        outFile2[barcode].write(line23)
                        outFile2[barcode].write(line24)
                else:
                    r1.readline()
                    r1.readline()
                    r1.readline()
                    r2.readline()
                    r2.readline()
                    r2.readline()
            else:
                r1.readline()
                r1.readline()
                r1.readline()
                r2.readline()
                r2.readline()
                r2.readline()
        line1 = r1.readline()
        line2 = r2.readline()
    for key in outFile:
        outFile[key].close()

def DualDemultiplex(run_id):
    fqdir = "/ufrc/cng/%s/fastq_files/rerun/"%(run_id)
    sampleSheet = "/ufrc/cng/%s/SampleSheet.csv"
    files = [f for f in os.listdir(fqdir) if f.endswith('_R1_001.fastq.gz')]
    for f in files:
        f2 = f[:-16] + "_R2_001.fastq.gz"
        cmd = "python ~/libHO/SeqFiles/script.py pipe.sortByBarcode %s %s %s"%(f,f2,fqdir)
        scriptOptions={'ppn':1,'jobname':'demultiplex','memory': '5gb','walltime':'36:00:00'}
        Mysbatch.launchJob(cmd,scriptOptions,verbose=True)"""
