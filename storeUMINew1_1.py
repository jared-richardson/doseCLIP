import gzip, sys
import operator

r1=open('/blue/berglund/jared.richardson/clip/Clip_Test3_2_1/trimmed/1_S1_R1_001.fastq', 'r')
r2=open('/blue/berglund/jared.richardson/clip/Clip_Test3_2_1/trimmed/1_S1_R2_001.fastq', 'r')
outfile=open('/blue/berglund/jared.richardson/clip/Clip_Test3_2_1/trimmed/1_S1_R1_001u.fastq', 'w')
outfile2=open('/blue/berglund/jared.richardson/clip/Clip_Test3_2_1/trimmed/1_S1_R2_001u.fastq', 'w')


line1 = r1.readline().strip("\n")
line12 = r1.readline().replace("\n", "")
line13 = r1.readline()
line14 = r1.readline().replace("\n", "")

line21 = r2.readline().strip("\n")
line22 = r2.readline().replace("\n", "")
line23 = r2.readline()
line24 = r2.readline().replace("\n", "")

d = {}
count1 = 0
count2 = 0

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

while len(line21) > 0:
	if line1[0] == "@" and line1[-8:] == "adapter1":
		umi = line12[-10:]
                umir = line22[:10]
                count1 += 1
                if umi not in d and (reverse_complement(umir) == umi):
                    d[umi] = []
                
                    line12u = line12[:-10]
                    line22u=line22[10:]
                    line14u= line14[:-10]
                    line24u=line24[10:]
                    outfile.write(line1 + "\n" + line12u + "\n" + line13 + line14u + "\n")
                    outfile2.write(line21 + "\n" + line22u + "\n" + line23 + line24u + "\n")
                    count2 += 1
        else:
             outfile.write(line1 + "\n" + line12 + "\n" + line13 + line14 + "\n")
             outfile2.write(line21 + "\n" + line22 + "\n" + line23 + line24 + "\n")
	line1 = r1.readline().strip("\n")
        line12 = r1.readline().replace("\n", "")
        line13 = r1.readline()
        line14 = r1.readline().replace("\n", "")
        line21 = r2.readline().strip("\n")
        line22 = r2.readline().replace("\n", "")
        line23 = r2.readline()
        line24 = r2.readline().replace("\n", "")

print count1
print count2	

