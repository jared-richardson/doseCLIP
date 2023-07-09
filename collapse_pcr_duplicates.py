import argparse, operator
import gzip, sys

def reverse_complement(kmer):
    """Takes a kmer of any length and returns the
        reverse complement of the nucleotides.

        Arguments:
        kmer -- Nucleotide sequence.

        Output:    
        reverse_complement -- Reverse complement of
            input sequence.
    """    
    alt_map = {'ins':'0'}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    for k, v in alt_map.items():
        # Lists out each nucleotide in the kmer.
        bases = list(kmer)
        # Grabs the complement of each base in the input kmer, puts in a
        # list, and then reverses the order of the list. Bases is a reversed
        # list iterator object.        
        bases = reversed([complement.get(base, base) for base in bases])
        # Joins nucleotides from iterator object. This is the reverse complement
        # of the input kmer.
        reverse_complement = ''.join(bases)
        return reverse_complement

def deduplicate_fastq(fastq_read1, fastq_read2,
                      fastq_read1_out, fastq_read2_out):
    """Iterates through FASTQ files and deduplicates or collapes
        PCR duplictes using the UMIs (randome 10mer) from the eCLIP paper. Outputs
        FASTQ files without PCR duplicates.

        Arguments:
        fastq_read1 -- First read of paired-end sequencing data 
            (FASTQ format).
        fastq_read2 -- Second read of paired-end sequencing data 
            (FASTQ format).

        Output:    
        fastq_read1_out -- First read output FASTQ file with PCR duplicates
            removed.
        fastq_read2_out -- Second read output FASTQ file with PCR duplicates
            removed.
    """
    # Initializes set used to store UMI sequences for each read. Uses set to
    # remove any reads with the same UMI sequence.
    umi_set = set()
    # Opens FASTQ files for processing,
    read1_opened = open(fastq_read1, 'r')
    read2_opened = open(fastq_read2, 'r')
    # Opens FASTQ files for output.
    read1_out = open(fastq_read1_out, 'w')
    read2_out = open(fastq_read2_out, 'w')
    # Reads all four lines from the two paired-end read input FASTQ Files.
    # These are used for iterated through all four lines in both files at once.
    line1 = read1_opened.readline().strip("\n")
    line12 = read1_opened.readline().strip("\n")
    line13 = read1_opened.readline().strip("\n")
    line14 = read1_opened.readline().strip("\n")
    line21 = read2_opened.readline().strip("\n")
    line22 = read2_opened.readline().strip("\n")
    line23 = read2_opened.readline().strip("\n")
    line24 = read2_opened.readline().strip("\n")
    # Loops through files until there are no more lines. Uses line21 but any of
    # the four lines would work. 
    while len(line21) > 0:
        # Checks to makes sure there is a proper heading for the line using the
        # "@" sign. Uses the "adapter1" keyword provided by Cutadapt to indicate
        # an adapter was removed and therefore the UMI sequence should be next.
        # Uses the ten nucleotides from both pair-end reads.
        if (line1[0] == "@") and (line1[-8:] == "adapter1"):
            umi_read1 = line12[-10:]
            umi_read2 = line22[:10]
            # Checks to make sure the second read also has a UMI that can confirm the read is
            # a PCR duplicate. If the second read does not confirm, the reads get written
            # out without removing UMI sequence and are not assumed duplicates.
            if (reverse_complement(umi_read2) == umi_read1):
                # Checks to see if UMI has been detected before. If it has, it does the
                # read gets discarded.
                if (umi_read1 not in umi_set):
                    # Adds UMI to set to check against later reads.
                    umi_set.add(umi_read1)
                    # Removes UMI sequences from reads. Removes nucleotide sequence and quality
                    # scores for these bases.
                    line12_clean = line12[:-10]
                    line22_clean = line22[10:]
                    line14_clean = line14[:-10]
                    line24_clean = line24[10:]
                    # Writes out modified read data.
                    read1_out.write(f"{line1}\n{line12_clean}\n{line13}\n{line14_clean}\n")
                    read2_out.write(f"{line21}\n{line22_clean}\n{line23}\n{line24_clean}\n")
            else:
                read1_out.write(f"{line1}\n{line12}\n{line13}\n{line14}\n")
                read2_out.write(f"{line21}\n{line22}\n{line23}\n{line24}\n")
        # If not adapter was trimmed from Cutadpat then the unmodified read set gets written
        # out. This is because of the inability to confirm the UMI sequence. 
        else:
            read1_out.write(f"{line1}\n{line12}\n{line13}\n{line14}\n")
            read2_out.write(f"{line21}\n{line22}\n{line23}\n{line24}\n")
        # Reads in next set of lines for the next paired-end reads.
        line1 = read1_opened.readline().strip("\n")
        line12 = read1_opened.readline().strip("\n")
        line13 = read1_opened.readline().strip("\n")
        line14 = read1_opened.readline().strip("\n")
        line21 = read2_opened.readline().strip("\n")
        line22 = read2_opened.readline().strip("\n")
        line23 = read2_opened.readline().strip("\n")
        line24 = read2_opened.readline().strip("\n")


def init_argparse():
    """Initiates the use of argparse. Returns parsed
        command line arguments.

        Arguments:
        See command line argument descriptions below.

        Output:
        Returns parser object. Object is used for function input.
    """
    # Script description.
    parser = argparse.ArgumentParser(description = "Collapses PCR duplicates using eCLIP \
                                                    UMI sequences for paired-end reads.")        
    # Command line arguments and descriptions.
    parser.add_argument("-r1", "--fastq_read1", action = "store", type = str, 
                        help="First read of paired-end sequencing data (FASTQ format).", 
                        required = True)
    parser.add_argument("-r2", "--fastq_read2", action = "store", type = str, 
                        help="Second read of paired-end sequencing data (FASTQ format).", 
                        required = True)  
    parser.add_argument("-r1o", "--fastq_read1_out", action = "store", type = str, 
                        help="First read output FASTQ file with PCR duplicates removed.", 
                        required = True)
    parser.add_argument("-r2o", "--fastq_read2_out", action = "store", type = str, 
                        help="Second read output FASTQ file with PCR duplicates removed.", 
                        required = True)                                                          
    return parser

def main():
    """Main script function. Executes other functions.

        Arguments:
        fastq_read1 -- First read of paired-end sequencing data 
            (FASTQ format).
        fastq_read2 -- Second read of paired-end sequencing data 
            (FASTQ format).

        Output:    
        fastq_read1_out -- First read output FASTQ file with PCR duplicates
            removed.
        fastq_read2_out -- Second read output FASTQ file with PCR duplicates
            removed.
    """
    # Uses argparse to organize command line arguments.
    args_to_parse = init_argparse()
    parsed = args_to_parse.parse_args()
    # Iterates through and executes main functions.
    deduplicate_fastq(parsed.fastq_read1, parsed.fastq_read2,
                      parsed.fastq_read1_out, parsed.fastq_read2_out)

# Executes main function.
if __name__ == "__main__":
    main()                       
