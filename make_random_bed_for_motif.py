import argparse
import random

def generate_regions(gtf_file, random_regions, replicates, output_file_prefix):
    """Generates random regions using input GTF file and input parameters.

    Arguemnts:
    gtf_file -- Gencode/Ensembl or other GTF file
        for the organism of interest. Should annotate
        the FASTA file used for alignment.
    random_regions -- Number of random regions to be
        generated and anlyzed for motifs.
    replicates -- Number of replicates to be generated
        for each random region set.
    output_file_prefix -- Prefix for output files.

    Output:
    bed_out -- BED file(s) of random generated regions.
        File name is output_file_prefix + "_random_#.bed"
    """
    # Generates list of exon events from GTF file.
    exon_list = get_exons_from_gtf(gtf_file)
    # Count used to count number of output files.
    file_count = 0 
    # List used to determine size of random region. 75 is weighted 4x
    # since it is the most common region size. 150 is weight 2x and
    # 300 is weighted 1x.
    region_size_list = [75, 75, 75, 75, 150, 150, 300]
    # First while loop iterates and generates a BED file for each
    # replicate set.
    while file_count < replicates:
        # Used to name output file.
        bed_out_pre = (output_file_prefix + "_random_" + str(file_count) + ".bed")
        bed_out = open(bed_out_pre, 'w') 
        file_count += 1
        # Used to count number of exons pulled from list. Randomly loops
        # through exon_list and grabs exons until the number of random
        # regions is reached.
        exon_count = 0
        while exon_count < random_regions:
            # Grabs random choice from list.
            exon_string = random.choice(exon_list)
            exon_sublist = exon_string.split("_")
            # Randomly chooses start or end position for random coordinates
            # to be based on.
            index_number = random.randrange(1, 2)
            random_cord = exon_sublist[index_number]
            # Randomly chooses a number of nucleotides to add or subtract.
            random_number_1 = random.randrange(-250, 250)
            # Randomly chooses a region size from region_size_list.
            region_size = random.choice(region_size_list)
            # Calculates start and end positions for random region.
            start = (int(random_cord) + random_number_1)
            end = (start + region_size)
            # Outputs region in BED format.
            bed_out.write(exon_sublist[0] + "\t" + str(start) + "\t" + str(end) +\
                          "\t" + "X" + "\t" + "100" + "\t" + exon_sublist[3] + "\n")
            exon_count += 1
                                   
def get_exons_from_gtf(gtf_file):
    """Generates list of exon events from input GTF file.

        Arguments:
        gtf_file -- Gencode/Ensembl or other GTF file
            for the organism of interest. Should annotate
            the FASTA file used for alignment.

        Output:
        exon_list -- List of exon events from GTF file.
            ["chrom_start_end_strand"]]
    """
    # Loops through GTF file and grabs exonic regions.       
    gtf_file_open = open(gtf_file, 'r')
    # List used to store exons. Exons will be pulled randomly
    # from this list.
    exon_list = []
    for line in gtf_file_open:
        if line[0] != "#":
            line_split = line.split("\t")
            # Filters for exonic regions anything with coordinates < 250.
            if (line_split[2] == "exon") and (int(line_split[3]) > 250):
                # Grabs chromosome, start, and end positions of exon.
                chrom, start, end, strand = (line_split[0], line_split[3], 
                                             line_split[4], line_split[6])
                exon_string = f"{chrom}_{start}_{end}_{strand}"
                exon_list.append(exon_string)
    gtf_file_open.close()
    return exon_list         

def init_argparse():
    """Initiates the use of argparse. Returns parsed
        command line arguments.

        Arguments:
        See command line argument descriptions below.

        Output:
        Returns parser object. Object is used for function input.
    """
    # Script description.
    parser = argparse.ArgumentParser(description = "Generates BED files of random regions to be used \
                                                    for motif analysis. Regions are used for normalization")        
    # Command line arguments and descriptions.
    parser.add_argument("-g", "--gtf_file", action = "store", type = str,
                        help="Gencode/Ensemlb or other GTF file for the organism of interest. \
                        Should annotate the FASTA file used for alignment.", 
                        required = True)
    parser.add_argument("-r", "--random_regions", action = "store", type = int,
                        help="Number of random regions to be generated and anlyzed for motifs.",
                        required = True)
    parser.add_argument("-rep", "--replicates", action = "store", type = int,
                        help="Number of replicates to be generated for each random region set.",
                        required = True)
    parser.add_argument("-o", "--output", action = "store", type = str,
                        help="Prefix for output files.",
                        required = True)                                                                 
    return parser

def main():
    """Main script function. Executes other functions.

        Arguments:
        gtf_file -- (-g) Gencode/Ensemlb or other GTF file
            for the organism of interest. Should annotate
            the FASTA file used for alignment.
        random regions -- (-r) Number of random regions to be
            generated and anlyzed for motifs.
        replicates -- (-rep) Number of replicates to be generated
            for each random region set.
        output_file_prefix -- (-o) Prefix for output files.             

        Output:    
        bed_out -- BED file(s) of random generated regions.
            File name is output_file_prefix + "_random_#.bed"

    """
    # Uses argparse to organize command line arguments.
    args_to_parse = init_argparse()
    parsed = args_to_parse.parse_args()
    # Iterates through and executes main functions.
    generate_regions(parsed.gtf_file, parsed.random_regions, parsed.replicates, parsed.output)

# Executes main function.
if __name__ == "__main__":
    main()