import argparse

def count_binding_regions(deseq2_files, region_bed_file, counts_file_csv):
    """Takes input DeSeq2 produced file and Piranha produced BED file,
        and finds the number of overlapping binding regions. Used for 
        finding overlapping binding regions for doseCLIP files and 
        another CLIP experimeint with no SM Input sample.

        Arguments:
        deseq2_files -- List of DeSeq2 produced differential expression
            file(s). Must be CSV file.
        region_bed_file -- BED file with binding regions. This file
            would be produced from a separate experiment where there
            where no SM Input samples.
        Output:    
        counts_file_csv -- Output counts file. Contains overlapping
            event counts for each input file.
    """
    # Opens outfile and writes out title.
    counts_file_csv_out = open(counts_file_csv, 'w')
    counts_file_csv_out.write("DeSeq2_filename,BED_filename,Overlapping_count\n")
    # Opens, loops through list of Deseq2 produced files, and iterates through, 
    # saving the chromosome coordinate variables for comparison.
    for deseq2_file in deseq2_files:
        # Sets count used for recording the number of overlapping events.
        count = 0
        with open(deseq2_file) as open_deseq2_file:
            file_name = deseq2_file.strip('"').split("/")[-1]
            for line in open_deseq2_file:
                line_sep = line.strip("\n").split(",")
                if line.find("~") != 0:
                    chrom_pre = line_sep[0].replace('"', '').split("~")
                    # Chromsome, event start coordinate 1, event start coordinate 2,
                    # strand.
                    chrom, cord1, cord2, strand = (chrom_pre[0], int(chrom_pre[1]), 
                                                   int(chrom_pre[2]), chrom_pre[3])
                    # Opens, loops through BED file, and iterates through, saving the 
                    # chromosome coordinate variables for comparison.    
                    with open(region_bed_file) as open_region_bed_file:
                        file_name2 = region_bed_file.strip('"').split("/")[-1]
                        for line2 in open_region_bed_file:
                            line_sep2 = line2.strip("\n").split("\t")
                            # Chromsome, event start coordinate 1, event start coordinate 2,
                            # strand.
                            chrom2, cord12, cord22, strand2 = (line_sep2[0], int(line_sep2[1]),
                                                              int(line_sep2[2]), line_sep2[5])
                            # Looks for any overlapping events from both files. If finds overlapping
                            # event, adds one to the count.
                            if (((chrom == chrom2) and (strand == strand2))
                                and (((cord1 <= cord12) and (cord12 <= cord2))
                                or ((cord1 <= cord22) and (cord22 <= cord2))
                                or ((cord1 <= cord12) and (cord22 <= cord2))
                                or ((cord1 >= cord12) and (cord22 >= cord2)))):
                                count += 1
        # Writes out file names and overlapping counts.                                                    
        counts_file_csv_out.write(f"{file_name},{file_name2},{count}\n")

def init_argparse():
    """Initiates the use of argparse. Returns parsed
        command line arguments.

        Arguments:
        See command line argument descriptions below.

        Output:
        Returns parser object. Object is used for function input.
    """
    # Script description.
    parser = argparse.ArgumentParser(description = "Takes input DeSeq2 produced file \
                                                    and Piranha produced BED file, \
                                                    and finds the number of overlapping \
                                                    binding regions. Used for finding \
                                                    overlapping binding regions for doseCLIP \
                                                    files and another CLIP experimeint with \
                                                    no SM Input sample.")        
    # Command line arguments and descriptions.
    parser.add_argument("-def", "--deseq2_file", action = "store", type = str, nargs='+', 
                        help = "DeSeq2 produced differential expression file(s). Must be CSV file",
                        required = True)
    parser.add_argument("-bed", "--region_bed_file", action = "store", type = str, 
                        help = "BED file with binding regions. This file would be produced \
                        from a separate experiment where there where no SM Input samples.",
                        required = True)
    parser.add_argument("-co", "--counts_file_csv", action = "store", type = str,  
                        help="Output counts file. Contains overlapping event counts for each \
                        input file.", 
                        required = True)                                                                              
    return parser

def main():
    """Main script function. Executes other functions.

        Arguments:
        deseq2_file -- List of DeSeq2 produced differential expression 
            file(s). Must be CSV file.
        region_bed_file -- BED file with binding regions. This file
            would be produced from a separate experiment where there
            where no SM Input samples.

        Output:    
        counts_file_csv -- Output counts file. Contains overlapping
            event counts for each input file.        
    """
    # Uses argparse to organize command line arguments.
    args_to_parse = init_argparse()
    parsed = args_to_parse.parse_args()
    # Iterates through and executes main functions.
    count_binding_regions(parsed.deseq2_file, parsed.region_bed_file,
                         parsed.counts_file_csv)

# Executes main function.
if __name__ == "__main__":
    main()