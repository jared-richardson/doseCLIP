import argparse

def add_sequence_output(bed_list, length, output):
    """Adds additional sequence based on user input to a BED file.

        Arguments:
        bed_list -- List of BED files of binding regions.
        length -- Length to add to both ends of all regions.
        output -- Output file directory.

        Output:
        Outputs the same BED files that were input with additional sequence
        added to the start and end of each region.    
    """
    # Iterates through each BED file in list.
    for bed in bed_list:
        # Opens BED file.
        with open(bed, "r") as bed_file:
            # Creates output file name.
            if output[-1] == "/":
                output_file_name = (output + bed.split("/")[-1].split(".")[0] + \
                                    "_additional_sequence.bed")
            else:
                output_file_name = (output + "/" + bed.split("/")[-1].split(".")[0] + \
                                    "_additional_sequence.bed")    
            # Opens output file.
            with open(output_file_name, "w") as output_file:
                # Iterates through each line in BED file.
                for line in bed_file:
                    # Splits line into list.
                    line_list = line.strip("\n").split("\t")
                    # Adds length to start and end of region. If the first
                    # cordinate is less than the length, the first cordinate
                    # is set to 1. Since the scaffold length is unknown, the
                    # second coordinate does not have this safeguard. Adjust
                    # the coordinate accordingly when converting to FASTA.
                    if int(line_list[1]) < length:
                        line_list[1] = "1"
                    else:    
                        line_list[1] = str(int(line_list[1]) - length)
                    line_list[2] = str(int(line_list[2]) + length)
                    # Writes line to output file.
                    output_file.write("\t".join(line_list) + "\n")


def init_argparse():
    """Initiates the use of argparse. Returns parsed
        command line arguments.

        Arguments:
        See command line argument descriptions below.

        Output:
        Returns parser object. Object is used for function input.
    """
    # Script description.
    parser = argparse.ArgumentParser(description = "Adds additional sequence based \
                                                    on user input to a BED file.")        
    # Command line arguments and descriptions.
    parser.add_argument("-b", "--bed", action = "store", type = str, nargs='+', 
                        help="BED file of binding regions.", required = True)
    parser.add_argument("-l", "--length", action = "store", type = int,
                        help="Length to add to both ends of all regions.",
                        required = True) 
    parser.add_argument("-o", "--output", action = "store", type = str,
                        help="Output file directory.", required = True)                                                                            
    return parser

def main():
    """Main script function. Executes other functions.

        Arguments:
        bed_list -- List of BED files of binding regions (-b/-bed).
        length -- Length to add to both ends of all regions (-l/-length).
        output -- Output file directory (-o/-output).

        Output:
        Outputs the same BED files that were input with additional sequence
        added to the start and end of each region.    

    """
    # Uses argparse to organize command line arguments.
    args_to_parse = init_argparse()
    parsed = args_to_parse.parse_args()
    # Iterates through and executes main functions.
    add_sequence_output(parsed.bed, parsed.length, parsed.output)

# Executes main function.
if __name__ == "__main__":
    main()