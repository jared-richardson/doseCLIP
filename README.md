# doseCLIP
This page contains the custom scripts required to process doseCLIP data. DoseCLIP is performed using the eCLIP protocol and uses their adapters/primers.  I will also 
include the commands for the use of open-source programs in the pipeline and the 
order of each steps.
## Step 1: Adapter trimming and PCR duplicate removal.
The first step is use Cutadapt to trim the reads. Please see the Cutadapt guide (https://cutadapt.readthedocs.io/en/stable/guide.html) for 
more information. The command I use for the first trim is  
as follows. The parameter -j should be set to the number of cores desired to be used when executing.
```
cutadapt -f fastq --match-read-wildcards --nextseq-trim=10 -q 10 --times 1  -e 0.05  -O 10 -m 18 -j 1 \
-y '{name}' \
-a adapter1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a adapter2=AGATCGGAAGAGCAC \
-g ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNCCTATAT -g ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNTGCTATT -g TCCGATCTNNNNNCCTATAT \
-g TCCGATCTNNNNNTGCTATT -g TTCCGATCTNNNNNCCTATA -g TTCCGATCTNNNNNTGCTAT -g CTTCCGATCTNNNNNCCTAT -g CTTCCGATCTNNNNNTGCTA \
-g TCTTCCGATCTNNNNNCCTA -g TCTTCCGATCTNNNNNTGCT -g CTCTTCCGATCTNNNNNCCT -g GCTCTTCCGATCTNNNNNCC -g CGCTCTTCCGATCTNNNNNC \
-g ACGCTCTTCCGATCTNNNNN -g CTCTTCCGATCTNNNNNTGC -g GCTCTTCCGATCTNNNNNTG -g CGCTCTTCCGATCTNNNNNT -g ACGCTCTTCCGATCTNNNNN \
-A ATATAGGNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -A AATAGCANNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-A ATATAGGNNNNNAGATCGGA -A AATAGCANNNNNAGATCGGA -A TATAGGNNNNNAGATCGGAA -A ATAGCANNNNNAGATCGGAA \
-A ATAGGNNNNNAGATCGGAAG -A TAGCANNNNNAGATCGGAAG -A TAGGNNNNNAGATCGGAAGA -A TCTTCCGATCTNNNNNTGCT \
-A AGGNNNNNAGATCGGAAGAG -A GGNNNNNAGATCGGAAGAGC -A GNNNNNAGATCGGAAGAGCG -A NNNNNAGATCGGAAGAGCGT \
-A GCANNNNNAGATCGGAAGAG -A CANNNNNAGATCGGAAGAGC -A ANNNNNAGATCGGAAGAGCG -A NNNNNAGATCGGAAGAGCGT \
-G GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNNNNNN -G GTGCTCTTCCGATCTNNNNNNNNNN \
-o output_fastq_file_read1.fastq \
-p output_fastq_file_read2.fastq \
input_fastq_file_read1.fastq \
input_fastq_file_read2.fastq \
> trim_metrics_file.txt
```
The next step is to remove the PCR duplicates. I developed a script to do this **collapse_pcr_duplicates.py**. The script also has a Pytest testing script 
(**test_collapse_pcr_duplicates.py**) and associated testing files. This can be used on your system by typing the `pytest` command in the script's directory (after 
installing pytest). The options for the **collapse_pcr_duplicates.py** script can be found using the command `python3 collapse_pcr_duplicates.py -h`. The script requires 
paired-end Illumina sequencing data.
**Step 2: Read alignment.
The next step is read alignment. Previously, I used to trim the reads a second time and then align. This used to be a requirement for many 
aligners, as the effectiveness of soft clipping varies by aligner. Minimap2 is highly effective at aligning regions that match genomic regions, 
making it so a second trimming is not necessary. The first trim would not be necessary either, if it not for the need to remove PCR duplicates. 
PCR duplicates are identified using the UMI sequence that itself is identified as being immediately adjacent to the adapter sequence. I used 
Minimap2 for aligning and the following is the command I use. The -t parameter is used to specify the number of threads desired to be used.
```  
minimap2 -t 1 -ax sr genomic_target_file.fasta input_fastq_file_read1.fastq input_fastq_file_read2.fastq > output_alignment_file.sam
```
For more information, see the Minimap2 manual at https://lh3.github.io/minimap2/minimap2.html.
## Step 3: Sorting and converting alignment files.
Before subsequent steps, the reads will need to be sorted. I used samtools to do this. The command is below. 
```samtools sort input_alignment_file.sam > output_alignment_file.sorted.sam
``` 
After this the SAM file will need to be converted to BAM format. I use the command below.
```samtools view -S -b input_alignment_file.sorted.sam > output_alignment_file.sorted.bam
```
For more information on samtools, please see the samtools manual- http://www.htslib.org/doc/samtools.html.
After using samtools, bedtools will need to be used to convert the BAM filse to BED files. The bedtools manual can be found here- 
https://bedtools.readthedocs.io/en/latest/. The command I use is below.
''' 
bedtools bamtobed -i output_alignment_file.sorted.bam > output_alignment_file.sorted.bed
'''
## Step 4: Using Piranha to determine read pileup locations and make custom annotation.
With the bed files, use Piranha to determine the bounderies of aligned reads. These regions will be used to build a custom GTF file that can be 
used to count read pileups in these regions. The Piranha manual can be found here- http://smithlabresearch.org/software/piranha/. The command I 
use for the program is below. The -z parameter is for the bin size and should be equal to the raw length of each input read.
```
Piranha -z 75 -u 0 -a 0.98 -s -o output_piranha_alignment_file.sorted.bed input_alignment_file.sorted.bed
```
Following using Piranha, a custom Python script is used to join all the regions of read pileups into a single annotation. Use the **join_binding_regions.py** script to do 
this. The script also has a Pytest testing script (**test_join_binding_regions.py**) and associated testing files. This can be used on your system by typing the `pytest` 
command in the script's directory (after installing pytest). The options for the **join_binding_regions.py** script can be found using the command `python3 
collapse_pcr_duplicates.py -h`. The script requires paired-end Illumina sequencing data. The script can take all the Piranha BED files and produce a single combined GTF 
file.
## Step 5: Counting the reads.
Now the reads are ready to be counted. To count the reads I use subread. The page for subread can be found here- 
https://subread.sourceforge.net/. The command I use is below. -T should be set to the number of threads desired to be used when running. All 
HITS-CLIP samples developed for each doseCLIP experiment should be counted together. The '-a' parameter should be the joined GTF file from the previous step. 
```
featureCounts -T 1 -p -O \
-a genomic_target_file.gtf \
-o output_counts_file.txt -t gene \
input_alignment_file1.sorted.bam \
input_alignment_file2.sorted.bam 
```
