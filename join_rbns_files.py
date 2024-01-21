import argparse

def get_file_name(file):
    """Gets file name from input argument.

        Arguments:
        file -- Input file string.

        Output:
        file_name -- File name without 
            directory structure and 
            file extension.
    """
    # Looks for "/" in directory structure. If found, splits
    # and grabs last element. Otherwise, just uses file name.
    if file.find("/") != -1:
        file_split = file.split("/")
        file_name = file_split[-1].split(".csv")[0]
    else:
        file_name = file.split(".csv")[0]
    return file_name          

def join_and_output_rbns(rbns_list, output_file):
    """Joins and outputs input motif data for RBNS
        style motif enrichment values.

        Arguments:
        rbns_list -- (-r) List of CSV files of RBNS style motif enrichment values.
        output_file -- (-o) Output directory.

        Output:
        CSV file of joined RBNS style motif enrichment values for 
            input files. First line contains file names used
            to generate motif enrichment values.
            File name- output_file + "all_joined_rbns_motif_enrichment.csv".
            Format - "Motif, Sample_Set1, Sample_Set2, ...". 
    """
    # Dictionary used to store motif enrichment values for all files. Stored
    # in value list in order of input files.
    # {"motif": [enrichment_value]}
    motif_dictionary = {}
    # List used to iterate through file_dictionary and grab files
    # and motifs in order.
    file_list = []
    # Gets sample number using length of rbns_list. Used to generate
    # 0 values for missing motifs.
    sample_number = len(rbns_list)
    # Opens each file and stores motif enrichment values in file_dictionary.
    for file in rbns_list:
        # Used to store file name.
        file_name = get_file_name(file)
        # Used to store file name in order.
        file_list.append(file_name)
        # Opens file for reading.
        open_file = open(file, 'r')
        # Loops through each line of file and stores motif enrichment
        # values in motif_dictionary.
        for line in open_file:
            if (line[0] != "#") and (line[0] != "M"):
                line_split = line.strip("\n").split(",")
                motif = line_split[0]
                enrichment_value = float(line_split[1])
                current_list = motif_dictionary.get(motif, [])
                current_list.append(enrichment_value)
                motif_dictionary[motif] = current_list
        open_file.close()
    # Adds "/" to end of output_file if not present.
    if output_file[-1] != "/":
        output_file = output_file + "/"
    # Opens output file for writing.
    output_file_name = (output_file + "all_joined_rbns_motif_enrichment.csv")
    output_file = open(output_file_name, 'w') 
    # Outputs title line.
    output_file.write("Motif")
    for file_name in file_list:
        output_file.write("," + file_name)
    output_file.write("\n")        
    # Iterate through motif_dictionary and outputs motif enrichment values.
    for motif in motif_dictionary:
        output_file.write(motif)
        # Count used to check against sample_number.
        count = 0
        for enrichment_value in motif_dictionary[motif]:
            count += 1
            output_file.write("," + str(enrichment_value))
        # Adds 0 values for missing motifs.
        if count < sample_number:
            for difference in range(sample_number - count):
                output_file.write(",0")    
        output_file.write("\n")
    output_file.close()

def init_argparse():
    """Initiates the use of argparse. Returns parsed
        command line arguments.

        Arguments:
        See command line argument descriptions below.

        Output:
        Returns parser object. Object is used for function input.
    """
    # Script description.
    parser = argparse.ArgumentParser(description = "Joins different sample concentrations of \
                                     final RBNS motif enrichment values into one file.")        
    # Command line arguments and descriptions.
    parser.add_argument("-r", "--rbns", action = "store", type = str, nargs='+', 
                        help="List of CSV file of RBNS style motif enrichment values.", 
                        required=True)      
    parser.add_argument("-o", "--output_file", action = "store", type = str,
                        help="Prefix for output files.",
                        required = True)                                                                          
    return parser

def main():
    """Main script function. Executes other functions.

        Arguments:
        rbns -- (-r) List of CSV files of RBNS style motif 
            enrichment values.
        output_file -- (-o) Output directory.

        Output:    
        CSV file of joined RBNS style motif enrichment values for 
            input files. First line contains file names used
            to generate motif enrichment values.
            File name- output_file + "all_joined_rbns_motif_enrichment.csv".
            Format - "Motif, Sample_Set1, Sample_Set2, ...". 
    """
    # Uses argparse to organize command line arguments.
    args_to_parse = init_argparse()
    parsed = args_to_parse.parse_args()
    # Iterates through and executes main functions.
    join_and_output_rbns(parsed.rbns, parsed.output_file)

# Executes main function.
if __name__ == "__main__":
    main()