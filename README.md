# doseCLIP
This page contains the custom scripts required to process doseCLIP data. DoseCLIP is performed using the eCLIP protocol and uses their adapters/primers.  I will also 
include the commands for the use of open-source programs in the pipeline and the 
order of each steps.
## Step 1: Read trimming and PCR duplicate demultiplexing.
The first step is use Cutadapt to trim the reads. Please see the Cutadapt manual for more information.  The command I use for the first trim is  as follows:
```
cutadapt -f fastq --match-read-wildcards --nextseq-trim=10 -q 10 --times 1  -e 0.05  -O 10 -m 18 \
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
-o /blue/berglund/jared.richardson/clip/clip_1/trimmed/1_S1_R1_001.fastq \
-p /blue/berglund/jared.richardson/clip/clip_1/trimmed/1_S1_R2_001.fastq \
/blue/berglund/jared.richardson/clip/clip_1/1_S1_R1_001.fastq \
/blue/berglund/jared.richardson/clip/clip_1/1_S1_R2_001.fastq \
> /blue/berglund/jared.richardson/clip/clip_1/TrimMetrics1
```
The next step is to remove the PCR duplicates. I developed a script to do this **collapse_pcr_duplicates.py**. The script also has a Pytest testing script 
(**test_pcr_duplicates**) and associated testing files. This can be used on your system by typing the `pytest` command in the script's directory. The options for the 
**collapse_pcr_duplicates.py** script can be found using the command `python3 collapse_pcr_duplicates.py -h`. The script requires paired-end Illumina sequencing data.
