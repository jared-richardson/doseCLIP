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
```
samtools sort input_alignment_file.sam > output_alignment_file.sorted.sam
``` 
After this the SAM file will need to be converted to BAM format. I use the command below.
```
samtools view -S -b input_alignment_file.sorted.sam > output_alignment_file.sorted.bam
```
For more information on samtools, please see the samtools manual- http://www.htslib.org/doc/samtools.html.
After using samtools, bedtools will need to be used to convert the BAM filse to BED files. The bedtools manual can be found here- 
https://bedtools.readthedocs.io/en/latest/. The command I use is below.
''' 
bedtools bamtobed -i output_alignment_file.sorted.bam > output_alignment_file.sorted.bed
'''
## Step 4: Using Piranha to determine read pileup locations and make custom annotation.
With the bed files, use Piranha to determine the bounderies of aligned reads. This only needs to be performed on the regular CLIP and the SM samples. 
These regions will be used to build a custom GTF file that can be used to count read pileups in these regions. The Piranha manual can be found 
here- http://smithlabresearch.org/software/piranha/. The command I use for the program is below. The -z parameter is for the bin size and should 
be equal to the raw length of each input read.
```
Piranha -z 75 -u 0 -a 0.98 -s -o output_piranha_alignment_file.sorted.bed input_alignment_file.sorted.bed
```
Following using Piranha, a custom Python script is used to join all the regions of read pileups into a single annotation. Use the **join_binding_regions.py** script to do 
this. The script also has a Pytest testing script (**test_join_binding_regions.py**) and associated testing files. This can be used on your system by typing the `pytest` 
command in the script's directory (after installing pytest). The options for the **join_binding_regions.py** script can be found using the command `python3 
collapse_pcr_duplicates.py -h`. The script requires paired-end Illumina sequencing data. The script can take all the Piranha BED files and produce a single combined GTF 
file. The script should be executed twice for each doseCLIP sample set. One GTF file should be generated using only the regular CLIP samples. This will be used for later
analyses involving normalized counts. A second GTF file should be generated using the regular CLIP files and the SM Input controls. This will be used to filter out binding
regions that are likely not authentically RBP bound regions. 
## Step 5: Counting the reads.
Now the reads are ready to be counted. Count all the CLIP and SM samples in two separate iterations. To count the reads I use subread. The page for subread can be found here- 
https://subread.sourceforge.net/. The command I use is below. -T should be set to the number of threads desired to be used when running. All 
HITS-CLIP samples developed for each doseCLIP experiment should be counted together. The '-a' parameter should be the joined GTF file from the previous step. One counts file
should be produced using the GTF that was generated using only the regular CLIP samples. This count file should only count the regular CLIP sample reads. A second count file
should be produced using the regular CLIP samples and the SM Input controls. The GTF file used for these counts should be the one that was generated using the SM Input 
controls and the regular CLIP samples. 
```
featureCounts -T 1 -p -O \
-a piranha_generated_gtf.gtf \
-o output_counts_file.txt -t gene \
input_alignment_file1.sorted.bam \
input_alignment_file2.sorted.bam 
```
## Step 6: Preparing the reads for DeSeq2.
After counting, the reads will need to be given titles for DeSeq2, formatted, and organized by sample type. I do this using Excel, but it can be done using a script, text 
editor, the command line, and/or a combination of the three. Subread outputs a few header lines with the rest of the file in txt format (tab deliminated). The minimum that 
needs to be done is removing the header lines, the extra columns removed, and the columns need titles. I also reformat to CSV but text or TSV can be input into R just as 
easily. Remove the header line(s) that indicates the processing commands, parameters, and samples. Insert a title row above the first line of counts. Remove the columns that indicate 
the chromosome (scaffold), the starting coordinate of the feature, the ending coordinate of the feature, and the strand. All that is needed is the gene name and the counts. 
Two counts files will be needed. One with just the normal CLIP samples and the other with the normal CLIP samples and the SM controls (these should have been generated in 
two separate processes in Step 5). The file with only the CLIP samples is used for determining counts for each binding region. The CLIP and SM contols counts file is used to filter out regions in the CLIP samples that do not show significantly more 
binding than observed in the SM samples. This removes false positives in the CLIP samples. The file should look something like this after (if using CSV format):
```
gene_id,sample_1,sample_2,sample_3,etc...
chr10~100283100~100283250~-g,3,4,5,etc...
```
Please note- the names used for samples will need to match the names in the phenotype file required for DeSeq2. The phenotype file is the files used by DeSeq2 to group 
samples together by type. The format I use for the file is as follows:
```
,condition,type
sample_1,x50,run_1
sample_2,x50,run_1
sample_3,x50,run_1
sample_4,x36,run_1
etc...
```
The condition is the sample type (in the case of doseCLIP it is protein concentration). Type is not required but can be useful for visualizations. These visualizations can 
help to understand if there were batch effects in the samples (i.e., sequencing certain samples together). Once the two counts files are prepared and the phenotype files 
prepared for each, the data can be input into DeSeq2 for differential expression analysis.
## Step 7: RNA-Seq data processing.
At this point it is advantageous to process the RNA-Seq data. Many of the steps are the same as with the CLIP data, except for the following differences. There of course are 
no SM controls for the RNA-Seq samples, so do not worry about these. The first step for the RNA-Seq samples (after ensuring they are of sufficient quality- this can be done 
with a tool like FastQC), is to align the reads. This can be performed using STAR for alignment. First an index of the target genomic sequence needs to be made. An example of 
the commands I use are below. The manual for STAR can be found here- https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf. The --runThreadN parameter should be 
set to the number of threads desired to be used and the --limitGenomeGenerateRAM parameter to the max RAM desired to be used.
```
STAR --runThreadN 32 --runMode genomeGenerate --genomeDir genomic_target_index_directory \
--sjdbGTFfile target_genomic_gtf_annotation.gtf \
--genomeFastaFiles genomic_target_file.fasta  \
--sjdbOverhang 74 \
--limitGenomeGenerateRAM=130000000000
```
After this, the RNA-Seq samples can be aligned. I use the following command.
```
STAR --runThreadN 32 --genomeDir genomic_target_index_directory \
--readFilesIn input_fastq_file_rnaseq_read1.fastq \
input_fastq_file_rnaseq_read2.fastq \
--outFileNamePrefix alignment_file_out.sorted.bam \
--outSAMtype BAM SortedByCoordinate
```
After performing the splicing analyses, subread can be used to get the gene counts from the BAM files. This can be performed using the same command as the CLIP/SM files, 
except for the GTF file used should be the one used for RNA-Seq alignment with STAR (should match the target genomic FASTA sequence). All RNA-Seq samples should be processed 
together. After this, the counts file should be prepared for DeSeq2, just like for the CLIP/SM samples. The RNA-Seq samples should have their own counts file, with only the 
RNA-Seq samples combined. The DeSeq2 phenotype file will also need to be prepared for the RNA-Seq samples. This file should contain all the RNA-Seq samples and again be 
separate. After preparing the files for DeSeq2, the splicing analysis can be performed. I use rMATS (turbo) for this. The rMATS manual can be found here- 
https://rnaseq-mats.sourceforge.net/rmats4.0.2/user_guide.htm. The command I use is below. The alternative splicing analysis should be performed comparing the uninduced 
samples versus each induced protein concentration. So for five induced sample sets and one uninduced sample set, there should be five rMATS analyses performed. The --nthread
parameter is how many threads are desired when running.
```
rmats.py --b1 x5_sample_set.csv \
--b2 uninduced_sample_set.csv \
--gtf target_genomic_gtf_annotation.gtf \
--od out_directory \
-t paired --readLength 76 --libType fr-firststrand --nthread 4
```
The x5_sample_set.csv should be a CSV file with the locations of the BAM files for each replicate in a sample set. For example, the contents of the x5_sample_set.csv file 
would be like below.
'''
x5_sample_1.sorted.bam,x5_sample_2.sorted.bam,x5_sample_3.sorted.bam
'''
The splicing analysis results can be saved for later analyses.
## Step 8: DeSeq2 analysis.
The **doseclip_deseq2_example.R** R script contains the needed R commands to be able to process the data. Edit the script for your samples and execute the commands. This 
should produce a variety of normalized counts files for downsteam analyses. There are also additional commands to produce plots like in the doseCLIP paper (not published yet). 
The main combined plot needs downstream analyses to be performed before it can be completed. After producing the required normalized counts files, proceed to the next steps.
## Step 9: Filtering the normalized counts 


