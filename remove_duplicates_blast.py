import argparse

def remove_duplicates(blast_out, output):
    """Removes duplicates using the first row region name
        in the input TSV file and then outputs non-duplicated
        TSV file.

        Arguments:
        blast_out -- List of BLAST produced output files.
        output -- Output directory.
        
        Output:
        blast_single_region -- Input BLAST files but with only one 
            region per organism comparison.
    """
    # Iterates through each BLAST output file.
    for file in blast_out:
        # Set used to check for duplicates.
        blast_single_set = set()
        # Opens BLAST output file.
        with open(file, "r") as blast:
            # Checks for "/" in output directory.
            if output[-1] != "/":
                output += "/"
            # Removes ".tsv" from file name.    
            if file.find(".tsv") != -1:
                file = file.replace(".tsv", "")
            # Creates output file name.
            output_file = (output + file.split("/")[-1].split(".")[0] + "_single_region.tsv")
            # Opens output file.
            with open(output_file, "w") as blast_single_region:
                # Iterates through each line in BLAST output file.
                for line in blast:
                    # Splits line by tab.
                    line = line.strip("\n").split("\t")
                    # If the first line of the file, write to output file.
                    if line[0] not in blast_single_set:
                        blast_single_region.write("\t".join(line) + "\n")
                        blast_single_set.add(line[0])
        
def init_argparse():
    """Initiates the use of argparse. Returns parsed
        command line arguments.

        Arguments:
        See command line argument descriptions below.

        Output:
        Returns parser object. Object is used for function input.
    """
    # Script description.
    parser = argparse.ArgumentParser(description = "Inputs final filtered DESeq2 files for all samples, \
                                                    removes lower concentration binding regions, and \
                                                    outputs normalized counts files for the input regions \
                                                    and the filtered regions.")        
    # Command line arguments and descriptions.
    parser.add_argument("-b", "--blast_out", action = "store", type = str, nargs='+',
                        help="List of BLAST produced output files", required = True)
    parser.add_argument("-o", "--output", action = "store", type = str, default = "",
                        help = "Output directory.", required = True)                                                                          
    return parser

def main():
    """Main script function. Executes other functions.

        Arguments:
        blast_out -- List of BLAST produced output files.
        output -- Output directory.  
            
        Output:    
        blast_single_region -- Input BLAST files but with only one 
            region per organism comparison.    
    """
    # Uses argparse to organize command line arguments.
    args_to_parse = init_argparse()
    parsed = args_to_parse.parse_args()
    # Iterates through and executes main functions.
    remove_duplicates(parsed.blast_out, parsed.output)

# Executes main function.
if __name__ == "__main__":
    main()