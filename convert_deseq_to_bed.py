import argparse

def convert_regions(deseq_files):
    """Converts DESeq2 output files into BED format

        Arguments:
        deseq_files -- List of DESeq2 produced files that need to be converted to BED format.
            Any DESeq2 file that contains gene names with the "~" separator should
            be able to be converted. Files will contain the same prefix as was input,
            and need to be in ".csv", ".txt", or ".tsv" format.

        Output:    
        BED files of the DESeq2 input are output.
    """
    # Iterates through each Deseq file, creates an output file, formats the event to BED,
    # and then outputs each event.
    for deseq_file in deseq_files:
        deseq_file_open = open(deseq_file, 'r')
        # Checks for normal file suffix. Places with ".bed". No if statements
        # as these just take time and do not have added benefit. If files contain
        # these suffixes other than the file suffix, this will change name other than
        # desired.
        deseq_file_out_pre = deseq_file.replace(".csv", ".bed")
        deseq_file_out_pre = deseq_file_out_pre.replace(".txt", ".bed")
        deseq_file_out_pre = deseq_file_out_pre.replace(".tsv", ".bed")
        deseq_file_out = open(deseq_file_out_pre, 'w')
        # Iterates through file and outputs in BED format.
        for line in deseq_file_open:
            if (line.find("~") != -1):
                if (line.find(",") != -1):
                    line_split = line.split(",")
                else:
                    line_split = line.split("\t")
                gene_split = line_split[0].split("~")
                # Need to grab only the strand which is the first character.
                gene_strand = gene_split[3][0]
                deseq_file_out.write(gene_split[0] + "\t" + gene_split[1] + "\t" + gene_split[2] +\
                                     "\t" + "X" + "\t" + "100" + "\t" + gene_strand + "\n")
        deseq_file_out.close()
        deseq_file_open.close()        
 
def init_argparse():
    """Initiates the use of argparse. Returns parsed
        command line arguments.

        Arguments:
        See command line argument descriptions below.

        Output:
        Returns parser object. Object is used for function input.
    """
    # Script description.
    parser = argparse.ArgumentParser(description = "Converts DESeq2 input files into BED format.")        
    # Command line arguments and descriptions.
    parser.add_argument("-d", "--deseq_file", action = "store", type = str, nargs='+', 
                        help="DESeq2 produced files that need to be converted to BED format. Files \
                              will be output to same directory",
                        required = True)                                                                              
    return parser

def main():
    """Main script function. Executes other functions.

        Arguments:
        deseq_files -- List of DESeq2 produced files that need to be converted to BED format.
            Any DESeq2 file that contains gene names with the "~" separator should
            be able to be converted. Files will contain the same prefix as was input,
            and need to be in ".csv", ".txt", or ".tsv" format.

        Output:    
        BED files of the DESeq2 input are output.
    """
    # Uses argparse to organize command line arguments.
    args_to_parse = init_argparse()
    parsed = args_to_parse.parse_args()
    # Iterates through and executes main functions.
    convert_regions(parsed.deseq_files)

# Executes main function.
if __name__ == "__main__":
    main()