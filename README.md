# doseCLIP Data Processing
This page contains the custom scripts required to process doseCLIP data (paper in preparation). DoseCLIP is performed using the eCLIP protocol and uses their adapters/primers. I will also include the commands for the use of open-source programs in the pipeline and the order of each step. Note- Python 3.6 or greater should be used to execute Python scripts.
## Step 1: Adapter trimming and PCR duplicate removal.
The first step is to use Cutadapt to trim the reads. Please see the Cutadapt guide (https://cutadapt.readthedocs.io/en/stable/guide.html) for 
more information. The command I use for the first trim is as follows. The parameter `-j` should be set to the number of cores desired to be used when executing.
```
cutadapt --match-read-wildcards --nextseq-trim=10 -q 10 --times 1  -e 0.05  -O 10 -m 18 -j 1 \
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
installing pytest using pip). The options for the **collapse_pcr_duplicates.py** script can be found using the command `python collapse_pcr_duplicates.py -h`. The script requires paired-end Illumina sequencing data. An example of how to run the script is below.
```
python /Users/jared.richardson/Desktop/doseclip/code/doseCLIP/collapse_pcr_duplicates.py \
-r1 1_S1_R1_001.fastq \
-r2 1_S1_R2_001.fastq \
-r1o 1_S1_R1_001_umi.fastq \
-r2o 1_S1_R2_001_umi.fastq
```
**Step 2: Read alignment.
The next step is read alignment. Previously, I used to trim the reads a second time and then align. This used to be a requirement for many 
aligners, as the effectiveness of soft clipping varies by aligner. Most modern aligners, including Subread, are highly effective at aligning regions that match genomic regions, 
making it so a second trimming is not necessary. The first trim would not be necessary either, if it not for the need to remove PCR duplicates. 
PCR duplicates are identified using the UMI sequence that itself is identified as being immediately adjacent to the adapter sequence. I used 
Subread for aligning and the following are the commands I use. The first command is to build the genome target index and the `-M` parameter is used to specify the amount of memory to use in megabytes. The second command is for read alignment. For this step, both the regular CLIP files and the SM controls should be aligned. The `-T` parameter is the number of threads to use. It is important to ensure that the PE orientation (`-S`) parameter is set correctly. For samples that were produced following the eCLIP protocol, it should match below (`-S fr`).
```
subread-buildindex -o index_directory -M ram_in_megabytes genomic_target_file.fasta
```
```
subread-align -t 1 -T threads_to_use -S fr -i index_directory \
-r input_fastq_file_read1.fastq \
-R input_fastq_file_read2.fastq \
-o output_alignment_file.bam
```
For more information, see the Subread website here https://subread.sourceforge.net/.
## Step 2: Sorting, converting, and filtering alignment files.
Before subsequent steps, the reads will need to be sorted. I used samtools to do this. The command is below. 
```
samtools sort input_alignment_file.bam > output_alignment_file.sorted.bam
```
For more information on samtools, please see the samtools manual- http://www.htslib.org/doc/samtools.html.
After using samtools, bedtools will need to be used to convert the BAM files to BED files. The bedtools manual can be found here- 
https://bedtools.readthedocs.io/en/latest/. The command I use is below.
``` 
bedtools bamtobed -i output_alignment_file.sorted.bam > output_alignment_file.sorted.bed
```
## Step 3: Using Piranha to determine read pileup locations and make custom annotation.
With the bed files, use Piranha to determine the bounderies of aligned reads. This only needs to be performed on the regular CLIP samples with the BED files. These regions will be used to build a custom GTF file that can be used to count read pileups in these regions. The Piranha manual can be found 
here- http://smithlabresearch.org/software/piranha/. The command I use for the program is below. The `-z` parameter is for the bin size and should 
be equal to the raw length of each input read.
```
Piranha -z 75 -u 0 -a 0.98 -s -o output_piranha_alignment_file.sorted.bed output_alignment_file.sorted.bed
```
Following using Piranha, a custom Python script is used to join all the regions of read pileups into a single annotation. Use the **join_binding_regions.py** script to do 
this. The script also has a Pytest testing script (**test_join_binding_regions.py**) and associated testing files. This can be used on your system by typing the `pytest` 
command in the script's directory (after installing pytest). The options for the **join_binding_regions.py** script can be found using the command 
`python join_binding_regions.py -h`. The script requires paired-end Illumina sequencing data. The script can take all the Piranha BED files and produce a single combined GTF 
file. The script should be executed once for each doseCLIP sample set. One GTF file should be generated using only the regular CLIP samples.
## Step 4: Counting the reads.
Now the reads are ready to be counted. Count all the CLIP and SM sample reads. To count the reads I use subread. The page for subread can be found here- 
https://subread.sourceforge.net/. The command I use is below. `-T` should be set to the number of threads desired to be used when running. All 
HITS-CLIP samples developed for each doseCLIP experiment should be counted together. The `-a` parameter should be the joined GTF file from the previous step. For me, it is easier to run the read count step in two different processes. One process that produces a counts file for only the regular CLIP files and the other with the regular CLIP files and the SM Input controls combined. The counts files will produce two different DESeq2 output file sets. The same GTF file should be used counting both sample sets. This should be the GTF that was generated from the previous step.
```
featureCounts -T 1 -p -O -s 1 \
-a piranha_generated_gtf.gtf \
-o output_counts_file.txt -t gene \
input_alignment_file1.sorted.bam \
input_alignment_file2.sorted.bam 
```
## Step 5: Preparing the reads for DESeq2.
After counting, the reads will need to be given titles for DESeq2, formatted, and organized by sample type. I do this using Excel, but it can be done using a script, text 
editor, the command line, and/or a combination of the three. Subread outputs a few header lines with the rest of the file in txt format (tab deliminated). The minimum that 
needs to be done is the header lines removed, the extra columns removed, and titles added to the columns. I also reformat to CSV but text or TSV can be input into R just as 
easily. Remove the header line(s) that indicates the processing commands, parameters, and samples. Insert a title row above the first line of counts. Remove the columns that 
indicate the chromosome (scaffold), the starting coordinate of the feature, the ending coordinate of the feature, and the strand. All that is needed is the gene name and the 
counts. Two counts files will be needed. One with just the normal CLIP samples and the other with the normal CLIP samples and the SM controls (these should have been 
generated in two separate processes in Step 5). The file with only the CLIP samples is used for determining counts for each binding region. The CLIP and SM contols counts 
file is used to filter out regions in the CLIP samples that do not show significantly more binding than observed in the SM samples. This removes false positives in the CLIP 
samples. The file should look something like this after (if using CSV format):
```
gene_id,sample_1,sample_2,sample_3,etc...
chr10~100283100~100283250~-g,3,4,5,etc...
```
Please note- the names used for samples will need to match the names in the phenotype file required for DESeq2. The phenotype file is the file used by DESeq2 to group 
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
help to understand if there were batch effects in the samples (i.e., sequencing certain samples together). Once the two counts files and the phenotype files 
prepared for each, the data can be input into DESeq2 for differential expression analysis.
## Step 6: RNA-Seq data processing.
At this point it is advantageous to process the RNA-Seq data. Many of the steps are the same as with the CLIP data, except for the following differences. There of course are 
no SM controls for the RNA-Seq samples, so do not worry about these. The first step for the RNA-Seq samples (after ensuring they are of sufficient quality- this can be done 
with a tool like FastQC- step not shown), is to align the reads. This can be performed using STAR for alignment. First an index of the target genomic sequence needs to be 
made. An example of the commands I use are below. The manual for STAR can be found here- https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf. The 
`--runThreadN` parameter should be set to the number of threads desired to be used and the --limitGenomeGenerateRAM parameter to the max RAM desired to be used.
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
together. After this, the counts file should be prepared for DESeq2, just like for the CLIP/SM samples. The RNA-Seq samples should have their own counts file, with only the 
RNA-Seq samples combined. The DESeq2 phenotype file will also need to be prepared for the RNA-Seq samples. This file should contain all the RNA-Seq samples and again be 
separate. After preparing the files for DESeq2, the splicing analysis can be performed. I use rMATS (turbo) for this. The rMATS manual can be found here- 
https://rnaseq-mats.sourceforge.net/rmats4.0.2/user_guide.htm. The command I use is below. The alternative splicing analysis should be performed comparing the uninduced 
samples versus each induced protein concentration. Therefore, for five induced sample sets and one uninduced sample set, there should be five rMATS analyses 
performed. The `--nthread` parameter is how many threads are desired when running.
```
rmats.py --b1 x5_sample_set.csv \
--b2 uninduced_sample_set.csv \
--gtf target_genomic_gtf_annotation.gtf \
--od out_directory \
-t paired --readLength 76 --libType fr-firststrand --nthread 4
```
The x5_sample_set.csv should be a CSV file with the locations of the BAM files for each replicate in a sample set. For example, the contents of the x5_sample_set.csv file 
would be like below.
```
x5_sample_1.sorted.bam,x5_sample_2.sorted.bam,x5_sample_3.sorted.bam
```
The splicing analysis results can be saved for later analyses.
## Step 7: DESeq2 analysis.
The **doseclip DESeq2_example.R** R script contains the needed R commands to be able to process the data. Edit the script for your samples and execute the commands. This 
should produce a variety of normalized counts files for downsteam analyses. There are also additional commands to produce plots like in the doseCLIP paper (not published yet). 
The main combined plot needs downstream analyses to be performed before it can be completed. After producing the required normalized counts files, proceed to the next steps.
## Step 8: Filtering the normalized counts 
After the normalized counts files for both the CLIP and SM Input dataset (all_clip_sm_normalized_counts.csv) and the regular CLIP dataset (all_clip_normalized_counts.csv) 
have been produced, as well as the CLIP and SM Input dataset SM filtered differential expression files (x50_clip_sm_filtered.csv, x36_clip_sm_filtered.csv, etc.) and the 
regular CLIP (in comparison to the uninduced) differential expression files (x50_vs_uni_significant.csv, x36_vs_uni_significant.csv, etc.), the regular CLIP produced files 
can be filtered. The binding regions also need to be annotated with the gene regions and sub-gene regions they are from. This can be done with the 
**filter_annotate_binding_regions.py** script. This script can be tested using the `pytest test_filter_annotate_binding_regions.py` command. For instructions on how to 
use the script, type `python filter_annotate_binding_regions.py -h`. Essentially, the script will annotate a CLIP and SM normalized counts file or differential 
expression file, if no regular CLIP file is also input. If the CLIP DESeq2 produced file is input with a regular CLIP file, then the CLIP and SM file will be annotated, the 
regular CLIP file will be annotated, and the regular file will be filtered using the SM filtered events. This will produce a SM filtered regular CLIP differential expression  
file but with the uninduced vs protein concentration differential expression information. The regular CLIP file normalized counts file will also need to be filtered and 
annotated. This can be done by concatenating all the different SM filtered protein concentrations (x50_clip_sm_filtered.csv, x36_clip_sm_filtered.csv, etc.) into a single 
file. Make sure to remove the title lines from the files before concatenating, except for the first title line. The filtered normalized counts are used for multiple later 
analyses. For a six protein concentration doseCLIP dataset (including uninduced), there should be twelve files produced. One with the regular CLIP normalized counts for all 
samples (`-o_cl` parameter in line 1 in the code example below). Five additional files that compared each protein concentration to the uninduced sample set for only the 
regular CLIP files (`-o_cl` paramter in line 2 in the code example below). Lastly, six SM/CLIP SM filtered events should be annotated (these require no DESeq2 produced 
regular CLIP files to be input when processing- `-o_sm` parameter in line 2 in the code example below). A GTF file is needed to annotate the binding regions and this 
should be the same GTF that was used for the RNA-Seq samples (should have been downloaded with the target genome FASTA file). Examples of how to execute thee 
**filter_annotate_binding_regions.py** script is below. There are two examples, for each of the file types just previously described.
```
python3 filter_annotate_binding_regions.py -sm all_clip_sm_filtered_concatenated.csv -gtf target_genomic_gtf_annotation.gtf -o_sm all_clip_sm_filtered_concatenated_annot.csv -clip 
all_clip_normalized_counts.csv -o_cl all_clip_normalized_counts_filt_annot.csv
python3 filter_annotate_binding_regions.py -sm x50_clip_sm_filtered.csv -gtf target_genomic_gtf_annotation.gtf -o_sm x50_clip_sm_filtered_annot.csv -clip x50_vs_uni_significant.csv -o_cl 
x50_vs_uni_significant_filt_annot.csv
``` 
In addition to the filtered and annotated output files from the filter/annotate step, counts files will be output for each input DESeq2 file. These counts files will have 
the ".annotation_counts.csv" suffix and will be output in the same directory as the other output files. These counts will be of the sub-gene features for each binding region.
These include- CDS, exon, intron, 5' UTR, 3' UTR, and intergenic. These counts are useful for determining the sub-gene regions the RBP is binding within the transcriptome.

To be able to make the MA plot with all datapoints, the five produced files need to be used. To use the code exactly like in the R script, the gene names will need to be 
copied and pasted into a new file with no column titles (or use a command like awk). These files are the ones that are imported in the R script. Alternatively, the original 
files can be imported into R and just the gene names extracted. You will need to write the code to do this, if desired.
## Step 9: Determining overlapping binding regions and annotating additional BED file(s).
In the doseCLIP paper (not published yet), MBNL1 binding regions are compared to a past Mbnl1 HITS-CLIP sample from a 2012 paper (Wang, et al., _Cell_). The HITS-CLIP 
sample was not made with an SM Input control making it so it could not be processed in the exact same way as the doseCLIP samples. The sample can have adaptors trimmed, PCR 
duplicates collapsed, reads aligned, and binding regions determined by Piranha. This ends up with a BED file containing aligned regions in the genome significantly above 
background (standard Piranha output). Although the sample will contain false positives, it is still useful for determining if the doseCLIP libraries were made successfully. 
The first step is to determine the overlapping regions between the past HITS-CLIP dataset and each filtered doseCLIP sample. This can be performed using the 
**count_binding_regions.py** script. See the script for execution instructions (`python count_binding_regions.py -h`). This script can be tested using the 
`pytest test_count_binding_regions.py` command. The previously annotated SM filtered DESeq2 produced files (the x50_clip_sm_filtered_annot.csv output file 
from the example above) should be checked for overlapping events. This means that for a doseCLIP sample set with six protein concentrations, including uninduced, the 
script should have six DESeq2 produced files input, with one BED file. The script will produce a text file with the counts of overlapping events for each DESeq2 produced 
file. These numbers can be used to make a figure showing the numbers of overlapping regions. Next, the Piranha produced BED file itself (from the past HITS-CLIP sample) can 
be annotated. This can be performed using the **filter_annotate_binding_regions_bed.py** script. See the script for execution instructions 
(`python filter_annotate_binding_regions_bed.py -h`). Again, this script can be tested using the `pytest test_filter_annotate_binding_regions_bed.py` command. The script 
needs a GTF file (should be the same one used to annotate the DESeq2 produced files in the previous step) and will produce an annotated BED file. This file can be used 
in the next step to determine the genomic distribution of the binding regions.
## Step 10 Determining Novel Regions Per Dose.
It was found in our studies that the characteristics of binding regions that are bound first (as determined after SM filtering and signifcant enrichment compared to the control, i.e., 'significant_filt_annot.csv' files), have different characteristics between lower and higher concentrations. These characteristics can be determined by filtering for Novel Regions Per Dose and comparing against the Comprehensive Regions (all filtered regions per concentration). To be able to filter these binding sites per concentration, the **remove_lower_regions.py** script can be used. See the script for execution instructions 
(`python remove_lower_regions.py -h`). Again, this script can be tested using the `pytest test_remove_lower_regions.py` command. The script needs all 'significant_filt_annot.csv' files, input in order (from lowest concentration to highest, NOTE- do not include the control or zero concentration!), the normalized counts file that contains the DESeq2 normalized counts for all samples (the 'all_clip_normalized_counts.csv' file, this should not include the SM samples!), a sample CSV file (contains the file to DESeq2 sample name), and a directory output location. The sample CSV file should be like the following format, with no title line. The sample names should be identical to what is location on the title line of the normalized counts file.
```
filename,deseq_sample_name1,deseq_sample_name2,deseq_sample_name3
filename,deseq_sample_name10,deseq_sample_name20,deseq_sample_name30
```
The script can be executed like the example below.
```
remove_lower_regions.py -de x5_vs_uni_significant_filt_annot.csv x10_vs_uni_significant_filt_annot.csv x20_vs_uni_significant_filt_annot.csv x36_vs_uni_significant_filt_annot.csv x50_vs_uni_significant_filt_annot.csv -n all_clip_normalized_counts.csv -csv sample_set.csv -o /homs/username/output/  
```
The script will produce two sets of outputs. One will contain the same output as the '_uni_significant_filt_annot.csv' files (the Comprehensive Binding Regions) but an additional last column will be added. This column will contain the average of the normalized counts. These files will have a '_count.csv' suffix. These files will need to be used for the later splicing analyses and can be used to make figures with normalized counts. The second set of files will contain a '_no_lower.csv' suffix and will contain the Novel Regions Per Dose for each sample. The Novel Regions Per Dose files are used for many downstream analyses. They will not be put in examples but will be noted when they need to used for processing by a script. Some steps, like the splicing steps, do no currently need the Novel Region Per Dose samples.
After this, the script **get_normalized_counts.py** should be executed using the '_uni_significant_filt_annot.csv' files created by **remove_lower_regions.py**. The script will can be executed like the example below.
```
python3 /Users/jared.richardson/Desktop/doseclip/code/doseCLIP/get_normalized_counts.py \
-de x5_vs_uni_significant_annotated_counts.csv \
10_vs_uni_significant_annotated_counts.csv \
x20_vs_uni_significant_annotated_counts.csv \
36_vs_uni_significant_annotated_counts.csv \
50_vs_uni_significant_annotated_counts.csv \
x5_vs_uni_significant_annotated_no_lower.csv \
x10_vs_uni_significant_annotated_no_lower.csv \
x20_vs_uni_significant_annotated_no_lower.csv \
x36_vs_uni_significant_annotated_no_lower.csv \
x50_vs_uni_significant_annotated_no_lower.csv \
-p data_set1 \
-o normalized_counts
```
The script will create one file with all average counts per each file combination and each single file ("_average_normalized_counts.csv" file suffix). In addition, files will be created containing each normalized count value for all singular files and files with counts based on the combinations of files. The files containing each normalized count value can be used for plotting all the values. If not all values are desired to be plotted, the averages can be used.
## Step 11: Determining enrichment of motifs (non-RNA Bind-N-Seq style).
First random BED files need to be generated that include regions around exonic regions for motif normalization. This can be performed with **make_random_bed_for_motif.py**. This script can be tested on your computer using `pytest test_make_random_bed_for_motif.py` and instructions for this script can be given with the following command- `python make_random_bed_for_motif.py -h`. The script takes a GTF file, an integer indicating the number of replicates to be generated, another integer indicating the number of regions per replicate to generate, and a file output prefix. I suggest generating at least three replicates using the average number of significant binding regions found for each dataset (perhaps around 3,000 regions). This will make it so the random distribution of your motifs of interest are likely generated in the random sets, enabling you to determine if your motifs of interest are enriched in your significant doseCLIP binding regions. The file prefix for these random BED files needs to include the string "random" in them. The motif analyzer script uses this string to determine which files are random and then uses these files for normalization. These random BED files will need to be converted to FASTA format like the filtered binding regions in the next step. The filtered binding regions, both verus the uninduced and SM Filtered can now be analyzed for motif enrichment. Motif enrichment informs us whether the RBP has nucleotide sequences that it has stronger binding interactions with. Before converting the files to BED files, it is advantageous at this point to make additional files that contain the 3' UTR binding regions for the 'uni_significant_filt_annot.csv' files and the the Novel Region Per Dose samples (i.e., '_no_lower.csv'). This can be performed using `grep` command with the phrase "3_UTR". The command can be performed like the example below. The 3 UTR files need to be converted to FASTA format but DO NOT need to be analyzed for motif information, at this point. They will be used for the binding region conservation step found later.
```
grep "3_UTR" x50_vs_uni_significant_filt_annot.csv > x50_vs_uni_significant_filt_annot_3utr.csv
```
The next step is to convert the DESeq2 produced binding events to BED format. This way Bedtools (website link is here: https://bedtools.readthedocs.io) can be used to grab the nucleotide sequences of these regions. Motif analysis can be performed on any genuine binding region, and therefore should be performed on the SM filtered binding regions (i.e., x50_clip_sm_filtered_annot.csv), the additionally filtered verus uninduced binding regions (i.e., x50_vs_uni_significant_filt_annot.csv), and the Novel Region Per Dose samples ('_no_lower.csv'). The first thing that needs to be done is the DESeq2 events need to be converted to BED format. This can be performed using the **convert_deseq_to_bed.py** script. This script can be tested on your system using Pytest as previosly mentioned. The script takes one or more DESeq2 produced files and outputs all the events in BED format in the file's directory. After this, Bedtools can be used to get the FASTA sequences for each binding region. To do this, the original FASTA sequence that was used for alignment will need to be used. This command is like the following:
```
bedtools getfasta -s -fi genomic_target_file.fasta -bed input.bed -fo output.fasta
```
After grabbing all the FASTA sequences, these sequences can be analyzed for motifs. This can be performed using the **random_enrichment_motif.py** script. Again, this script can be tested on your computer using `pytest test_random_enrichment_motif.py` and instructions for this script can be given with the following command- `python random_enrichment_motif.py -h`. The file takes a list of FASTA files that need to be analyzed for motif enrichment. This list should include the random files preiously generated. These random files will be used for normalization (these random files need "random" in the name). The python script also takes a kmer size. If looking for MBNL1 binding regions, this value will be four. There should also be an output file prefix and if looking for motifs that are not for MBNL1, then a list of motifs can be input, as well as a title for this set of motifs. These last two parameters are optional and do not need to be input, if looking for YGCYs for MBNL1. An example command for the script for MBNL1 can be found below. All FASTA files for each doseCLIP set can be input together. This step should also be performed on the Novel Region Per Dose files!
 ```
random_enrichment_motif.py -fl x50_vs_uni_significant_filt_annot.fasta x50_clip_sm_filtered_annot.fasta random_1.fasta random_2.fasta random_3.fasta -k 4 -o mbnl1_doseclip  
```
The output includes three files. One file that has the motif enrichment per 100 nucleotides for each binding region in the dataset ("_motifs_per_region.csv") suffix, a file that shows the normalized enrichment for all motifs ("all_motifs.csv" suffix), and then a file that includes the normalized enrichment for only the motifs of interest ("total_motifs.csv" suffix). Please see the Python script for more details on the output files. Next, another script can be executed to output motif secondary features, other than the primary motif enrichment data. This script takes the same input as **random_enrichment_motif.py** (separated because of the length of the code and separate functionality), except for it does not need the random input FASTA files. The script is **motif_secondary_features.py** and can be tested or help given like above. An example of the script is below. This step should also be performed on the Novel Region Per Dose files!
 ```
motif_secondary_features.py -fl x50_vs_uni_significant_filt_annot.fasta x50_clip_sm_filtered_annot.fasta ... -k 4 -o mbnl1_doseclip 
```
The script will produce four files. Two files are produced for each type of output since there are so many columns per output. One set of output contains the secondary motif statistics per file (suffix- '_secondary_motif_stats_per_file_1.csv') and per sample but on a per region basis (suffix- _secondary_motif_stats_per_region_1.csv). The information from these files can be directly used to make final figures for the secondary motif statistics.
## Step 12: Binding region and splicing analysis.
After the motifs have been analyzed, the splicing events (excluding MXE events, the script will need to be modified for these) can be related to the binding regions. This can be performed using the **match_binding_to_splicing.py** script. This script will also produce a filtered splicing file that contains events that have an FDR of 0.10 or less and a delta PSI of 0.10 or greater. The script can be tested using `pytest test_match_binding_to_splicing.py` and instructions for this script can be given with the following command- `python match_binding_to_splicing.py -h`. The script needs the RMATS produced "SE.MATS.JCEC.txt" file (for SE splicing event analysis), an annotated and filtered DESeq2 file (i.e., "x50_vs_uni_significant_filt_annot.csv" or "x50_clip_sm_filtered_annot.csv"), a motifs per region CSV file ("_motifs_per_region.csv") from the previous step, the distance from the binding region endpoint to the regulated exon, and a prefix for the output file (should contain a separator like an underscore at the end). The files should be all be for the same sample (i.e., 50x). An example command for the script can be found below. None of the splicing steps currently need the Novel Region Per Dose samples.
```
match_binding_to_splicing.py -r SE.MATS.JCEC.txt -b x50_vs_uni_significant_filt_annot.csv -m x50_vs_uni_significant_filt_annot_motifs_per_region.csv -d 1000 -o /home/username/output/x50_vs_uni_significant_filt_annot_se_ 
```
The script will produce the filtered RMATS splicing file (suffix- "filtered_SE.MATS.JCEC.txt") and a file (suffix- "splicing_binding.csv") containing information about the binding regions and splicing events within the distance specified. These files can be making splicing figures and the next step. The next step is to use **analyze_splicing_sample_set.py** script to find all filtered splicing events with proximal binding regions, all proximal, common sample splicing events, and splicing events that had more than one discrete binding region within the distance specified in **match_binding_to_splicing.py**. The script can be tested using `pytest test_analyze_splicing_sample_set.py` and instructions for this script can be given with the following command- `python analyze_splicing_sample_set.py -h`. The script needs the output (the 'filtered_SE.MATS.JCEC.txt' and the 'splicing_binding.csv' files for all samples) from the previous script to run. In addition, the sample_csv_file  used for the `remove_lower_regions.py` script needs to be input, as does the unfiltered  "SE.MATS.JCEC.txt" files from the previous step, and the DESeq2 created normalized counts file, containing all the raw normalized counts. All sample files except for the control or zero concentration (this should not have a splicing file anyways since it was used as the control) should be input for the entire doseCLIP sample set. Also, all samples should be in the same order in each file list and this should be the order in the sample_csv_file. The script will produce two files and another set of files. One file that gives the splicing events that are commonly signifcant in all samples  ('_all_filtered_SE.MATS.JCEC.txt') and a file that gives the splicing events and binding region information for splicing events where multiple binding regions were found '_multiple_binding_splicing.csv'. The set of files will have the suffix '_multiple_file_regions.csv' and one output will be created for each file. The output will contain binding and splicing information for all splicing events that had binding regions within the specified distance and were found in multiple sample sets. Analyzing events found in multiple samples has shown to produce more robust results. The script can be executed like the example below>
 ```
analyze_splicing_sample_set.py -s x5_vs_uni_significant_filt_annot_se_filtered_SE.MATS.JCEC.txt x10_vs_uni_significant_filt_annot_se_filtered_SE.MATS.JCEC.txt ... -b x5_vs_uni_significant_filt_annot_se_splicing_binding.csv x10_vs_uni_significant_filt_annot_se_splicing_binding.csv ... -o /home/username/output/all_splicing
```
This should be the last script needed to be executed before making figures related to splicing. All the splicing script outputs can be easilly processed to make figures.
## Step 13: Binding Region Conservation Analysis
Using the FASTA files from Step 11, the regular SM Filtered files (i.e., x50_clip_sm_filtered_annot.fasta), execute the regions through BLAST to determine the conserved regions with specific organisms. For the original doseCLIP study, the experiment was performed in a mouse derived cell line and conservation was analyzed for rat, Rhesus monkey, human, cow, and zebrafish. The taxonomic IDs for each of the organisms as well as their scientific names are listed below.
```
Scientific Name, Formal Name, Taxonomic ID
Rattus norvegicus, rat, 10116
Macaca mulatta, Rhesus monkey, 9544
Homo sapiens, human, 9606
Bos taurus, cow, 9913
Danio rerio, zebrafish, 7955
```
Depending on the organism studied, different organisms might need to be analyzed. Zebrafish were analyzed for doseCLIP because they are one of the "lowest level" metazoans that has all three paralogs of _mbnl_. The conservation did not end up being signifcant, so no figures included zebrafish. There are three ways to use BLAST, use the graphic interface, using the elasticblast resource to submit jobs, or using a local version of BLAST (recommended). If you want to use a local version of BLAST and do not have it currently, instructions for how to set it up can be found here- https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata. If using a local version, it is advantageous to write a loop or add all the commands to a bash script. I wrote an example for this in the **blast_submission.sh** script. The script requires a directory filled with ".fasta" files to analyze and a taxonomic ID. The script will output the files in the input directory and they will contain an "_blast.tsv" extension. An example of how to run the script against rat sequence data is below. See the script for more details.
```
./blast_submission.sh input_directory 10116
```
Alternatively, the BLAST elasticblast can be used to execute searches. The elasticblast instructions can be found here- https://blast.ncbi.nlm.nih.gov/doc/elastic-blast/elasticblast.html. I believe the elasticblast option requires spending for compute through a Cloud provider, so keep that in mind. The most accessible way to run BLAST for many users will be using the graphic interface. To do this go to the BLAST website- https://blast.ncbi.nlm.nih.gov/Blast.cgi. Click "Nucleotide BLAST". Enter a file for a specific sample under "Or, upload file". Under the "Organism" secion, enter the taxomic ID for the organism to search against. The full title for the organism should appear in the pull down menu. Click the full name. Finally, under the "Optimize for" section, click "Somewhat similar sequences (blastn)".  After this you can submit the search. The jobs can take a while to execute. It is best not to submit too many jobs at once as NCBI monitors this and will cancel submissions if too many are submitted at the same time (I am not sure what the current limits are). After this use **remove_duplicates_blast.py** with the output files to remove duplicates in each BLAST produced file. The script can be tested like the other Python scripts and help can be given with the help parameter. The script only needs one or more BLAST files as input and an output directory. The script also requires the annotated SM Filtered DeSeq2 output CSV files produced in Step 8. They are used to count the number of each type of binding regions in the output file. They should be in the same order as the BLAST input files.  The script can be executed like the example below.
```
remove_duplicates_blast.py -b 50x_blast_output.tsv -d x50_clip_sm_filtered_annotated.csv  -o output_directory
```
A counts file will be produced from this script ("_blast_counts_file.tsv"), as well as the the conserved, non-duplicate binding regions. The files produced from this script can be used directory for making the figures shown in the doseCLIP paper.
## Step 14: Structure Analysis
For the structure analysis ViennaRNA RNAfold will need to be used. The instructions for how to setup the program can be found here- https://www.tbi.univie.ac.at/RNA/. Before running RNAfold, the BED files used to create the FASTA files in Step 11 will need to be modified to increase the length of each region. The length is increased to get more of the actual RNA secondary structure features, as the regions themselves likely form secondary structures with the nucleotides surrounding them. I suggest adding at least 500 nucleotides to each side of the binding regions. Ideally, you would use the entire RNA transcript, but computationally this gets more difficult to process and the software is likely not able to predict such large secondary structures accurately. I wrote a script called **add_additional_sequence.py** to do this. The script can be tested like the examples before and further instructions for executing can be given by using the `-h` argument with the script. The script requires a list of BED files where processing is needed, and integer length of the number of nucleotides to be added to both ends of the region, and an output file directory. An example of how to run the script can be found below.
```
python add_additional_sequence.py -b x50_vs_uni_significant_filt_annot.bed x36_vs_uni_significant_filt_annot.bed -l 500 -o output_directory
```
All steps in Step 14 should be exectuted with the additionally filtered verus uninduced binding regions (i.e., x50_vs_uni_significant_filt_annot.fasta) and the Novel Region Per Dose samples ('_no_lower.fasta'.) After this, the output BED files should be converted to FASTA files. This can be performed with Bedtools like the example below and in Step 11.
```
bedtools getfasta -s -fi genomic_target_file.fasta -bed input.bed -fo output.fasta
```
After this ViennaRNA RNAfold should be ran on the FASTA files. I wrote a script that will perform this step on all FASTA files in a directory. Please edit the location of the ViennaRNA executables if alias does not match the script (and/or is installed system wide).
```
./get_rna_structure.sh fasta_directory
```
The output ("_structure_fasta") contains the FASTA information with an additional line that indicates whether the nucleotides in the position are likely to be paired with another nucleotide in the sequence or not. Using the **motif_secondary_features_structure.py** script, these files can be used to produce counts of the motifs and their surrounding sequence to determine the percentages of paired bases. The file can be tested and help given as previously described with Python scripts. An example of how to run the script is below. The script needs the same input as **motif_secondary_features.py** but only outputs two files. One file containing the structural information per region "secondary_structure_per_region.csv" and the other containing the information per file ("secondary_structure_per_file.csv"). These two files can be used directly to make the final figures.
```
python motif_secondary_features_structure.py -fl x50_vs_uni_significant_filt_annot_structure_fasta -k 4 -o output_directory
```
## Step 14: RNA Bind-N-Seq Style Motif Enrichment Analysis
The data can also be analyzed in the enrichment style of RNA Bind-N-Seq (RBNS), from the 2014 (Lambert et al., _Mol Cell_) paper. This type of analysis has never been performed on seuuencing libraries that were not _in vitro_, so it likely enrichment values are directly comparable, but are still informative to relate in some aspects. The RBNS enrichment values can be calculated using the **rbns_motif_enrichment.py** script. The script can be tested and help given as mentioned previously. The script needs the trimmed adn UMI removed FASTQ file (the final FASTQ files from Step 1) for each sample and its SM Input control. The script also needs an output directory location, a sample prefix, and kmer size (this should be seven for MBNL, as this is the kmer size used in the RBNS paper). The sample sets and their respective SM Inputs should be input all together and in order (i.e., sample1, sample2, sample3 - sm_sample1, sm_sample2, sm_sample3). The script will need to be executed once per sample set (so if six different concentrations, including uninduced, it should be executed six times). Uninduced samples should be analyzed as well. The script can be executed like the example below.
```
python rbns_motif_enrichment.py -f sample1.fastq sample12.fastq sample13.fastq -fs sm_sample1.fastq sm_sample12.fastq sm_sample13.fastq -k 7 -o output_directory -p uninduced_samples
```
After running the script on all samples, **join_rbns_files.py** can be executed with all the output files so join all the data into a single file for easy data retrieval. The script can be tested and help given as mentioned previously. All sample concentration RBNS motif files should be input together. The script can be executed like the example below.
 ```
python join_rbns_files.py -r uninduced_rbns_motif_enrichment.csv x5_rbns_motif_enrichment.csv x10_rbns_motif_enrichment.csv x20_rbns_motif_enrichment.csv x36_rbns_motif_enrichment.csv x50_rbns_motif_enrichment.csv -o output_directory
```
This is currently all the steps needed to recreate the analyzes performed in the original doseCLIP paper (not published yet). Further code will be added so that the pipeline is a single execution, with only a configuration file needed to be filled to complete. This will be a separate repository. In additon, motif analysis steps using machine learning will also be added. Currently, I am unsure if I will add the steps here or just point to the open-source software used for the analyses.