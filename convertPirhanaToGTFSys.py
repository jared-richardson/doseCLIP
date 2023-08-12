import sys

#r1=open('/blue/berglund/jared.richardson/clip/clipTrial9_1/bam/s2piranha', 'r')
outfile=open(sys.argv[2], 'w')
count = 0
d = {}

with open(sys.argv[1]) as f:
        for x in f:
            x1 = x.strip("\n")
            cols2=x1.split('\t')  
            #if cols2[0] + "\t" + cols2[1] + "\t" + cols2[2] not in d:
            outfile.write(cols2[0] + "\t" + "pir" +"\t"+ "gene" +"\t" +cols2[1] + "\t" + cols2[2] + "\t" +cols2[4] + "\t" +cols2[5]+"\t"+cols2[6] \
                         +"\t"+'gene_id "' +cols2[0] + "~"+cols2[1] + "~"+ cols2[2]+"~"+cols2[5] +'g"; ' + 'gene_version "1"; ' + "\n")
            outfile.write(cols2[0] + "\t" + "pir" +"\t"+ "transcript" +"\t" +cols2[1] + "\t" + cols2[2] + "\t" +cols2[4] + "\t" +cols2[5]+"\t"+cols2[6] \
                         +"\t"+'gene_id "' +cols2[0] + "~"+cols2[1] + "~"+ cols2[2]+"~"+cols2[5] + '"; ' + 'gene_version "1"; ' + \
                         'transcript_id "' +cols2[0] + "~"+cols2[1] + "~"+ cols2[2]+"~"+cols2[5] + 't"; ' + "\n")
            #d[cols2[0] + "\t" + cols2[1] + "\t" + cols2[2]] = []
            #print(cols2)
