import argparse

def find_output_all_splicing(filtered_splicing, output):
    """Uses file list to find splicing events present in all
        files and outputs them in the same format.

        Arguments:
        filtered_splicing -- rMATS style output file that is
            filtered by delta PSI and an FDR of 10% or less.
            Contains the 'filtered_SE.MATS.JCEC.txt' suffix.
        output -- Prefix for output file.

        Output:
        splicing_events_all_files_out -- Filtered rMATS style 
            output file with splicing events from all input files.
            The orginal filename is added to the last column. 
            File name- prefix + "_all_filtered_SE.MATS.JCEC.txt".
    """
    # Dictionary key is the splicing event, and the value is the list 
    # of lines.
    splicing_events_all_files = {}
    # List used to store file names.
    file_names = []
    # Used to save the title line for output.
    title_line = ""
    # Opens output file for writing.
    splicing_events_all_files_out = open(output + "_all_filtered_SE.MATS.JCEC.txt", "w")
    # Loops through files, and adds splicing events to a dictionary. If splicing
    # event is in all files, outputs the splicing event and the filename.
    for file in filtered_splicing:
        # Gets filename. If directory name is included, removes it.
        if file.find("/") != -1:
            file_name = file.split("/")[-1]
        else:
            file_name = file
        file_names.append(file_name)    
        with open(file, "r") as f:
            for line in f:
                line_clean = line.strip("\n")
                if line_clean[0] != "I":
                    line_list = line.split("\t")
                    # splicing_event = (chromosome, exon-start, exon-end, up-start, up-end, down-start, down-end, strand)
                    splicing_event = ("%s:%s:%s:%s:%s:%s:%s:%s" % (line_list[3], line_list[5], line_list[6], 
                                                                   line_list[7], line_list[8], line_list[9], 
                                                                   line_list[10], line_list[4]))
                    # If splicing event is not in dictionary, adds it.
                    if splicing_event not in splicing_events_all_files:
                        splicing_events_all_files[splicing_event] = [line_clean]
                    # If splicing event is in dictionary, adds line to list.
                    else:
                        splicing_events_all_files[splicing_event].append(line_clean)
                else:
                    title_line = line_clean                     
    # Outputs title line to file.
    splicing_events_all_files_out.write("%s\tFile_Name\n" % (title_line))
    # Checks for splicing events by using the length of the list. If it matches
    # the number of files, it is in all files. These will be written to the output file.
    for splicing_event in splicing_events_all_files:
        # Count used to output correct filename.
        file_count = 0 
        if len(splicing_events_all_files[splicing_event]) == len(filtered_splicing):
            for line in splicing_events_all_files[splicing_event]:
                splicing_events_all_files_out.write("%s\t%s\n" % (line, file_names[file_count]))
                file_count += 1

def find_output_multiple_binding(splicing_binding, output):
    """Uses file list to find splicing events with multiple
        binding regions around a single spliced exon.

        Arguments:
        splicing_binding -- CSV file with splicing and binding
            information for each region. Produced by the 
            match_binding_to_splicing.py script.
        output -- Prefix for output file.

        Output:
        multiple_binding_splicing_out -- CSV file with splicing and
            multiple binding regions around a single spliced exon.
            The distance of the double binding regions was specified
            in the match_binding_to_splicing.py script. 
            File name- prefix + "_" + filename + "_multiple_binding_splicing.csv".
    """
    # Dictionary key is the splicing event, and the value is the list 
    # of lines.
    splicing_binding_all_files = {}
    # List used to store file names.
    file_names = []
    # Used to save the title line for output.
    title_line = ""
    # Opens output file for writing.
    multiple_binding_splicing_out = open(output + "_multiple_binding_splicing.csv", "w")
    # Loops through files, and adds binding/splicing events to a dictionary. If more than
    # one binding region is found, outputs the splicing event and binding regions.
    for file in splicing_binding:
       # Gets filename. If directory name is included, removes it.
        if file.find("/") != -1:
            file_name = file.split("/")[-1]
        else:
            file_name = file
        file_names.append(file_name)
        with open(file, "r") as open_file:
            for line in open_file:
                line_clean = line.strip("\n")
                if line_clean[0] != "R":
                    line_list = line.split(",")
                    # splicing_event = (chromosome, exon-start, exon-end, strand, 
                    # delta-PSI)
                    # NOTE- Since the splicing event does not contain the up-stream 
                    # region and the down-stream region, it could result in a false 
                    # positive. This is unlikely, but should be noted. This will likely 
                    # be changed in future revisions (revisions need to be performed 
                    # with match_splicing_to_binding.py). 
                    splicing_event = ("%s:%s:%s:%s:%s" % (line_list[2], line_list[3], line_list[4], 
                                                          line_list[5], line_list[7]))
                    # If splicing event is not in dictionary, adds it.
                    if splicing_event not in splicing_binding_all_files:
                        splicing_binding_all_files[splicing_event] = [line_clean + "," + file_name]
                    # If splicing event is in dictionary, adds line to list.
                    else:
                        splicing_binding_all_files[splicing_event].append(line_clean + "," + file_name)
                else:
                    title_line = line_clean
    # Writes title line to file.
    multiple_binding_splicing_out.write("%s,Filename\n" % (title_line))                      
    # Outputs multiple binding regions to file.
    for splicing_event in splicing_binding_all_files:
        if len(splicing_binding_all_files[splicing_event]) > 1:
            for line in splicing_binding_all_files[splicing_event]:
                multiple_binding_splicing_out.write("%s\n" % (line)) 


def init_argparse():
    """Initiates the use of argparse. Returns parsed
        command line arguments.

        Arguments:
        See command line argument descriptions below.

        Output:
        Returns parser object. Object is used for function input.
    """
    # Script description.
    parser = argparse.ArgumentParser(description = "Uses files from match_binding_to_splicing.py \
                                                    to output splicing events that are in all files \
                                                    and multiple binding region splicing events.")
    # Command line arguments and descriptions.
    parser.add_argument("-s", "--filtered_splicing", action = "store", type = str, nargs='+',
                        help="rMATS style output file that is \
                              filtered by delta PSI and an FDR of 10 percent or less. \
                              Contains the 'filtered_SE.MATS.JCEC.txt' suffix.",                          
                        required = True)
    parser.add_argument("-b", "--splicing_binding", action = "store", type = str, nargs='+',
                        help="CSV file with splicing and binding \
                              information for each region. Produced by the \
                              match_binding_to_splicing.py script.",
                        required = True)
    parser.add_argument("-o", "--output", action = "store", type = str,
                        help="Prefix for output files.",
                        required = True)                                                                 
    return parser

def main():
    """Main script function. Executes other functions.

        Arguments:
        filtered_splicing -- rMATS style output file that is
            filtered by delta PSI and an FDR of 10% or less.
            Contains the 'filtered_SE.MATS.JCEC.txt' suffix.
        splicing_binding -- CSV file with splicing and binding
            information for each region. Produced by the 
            match_binding_to_splicing.py script.      
        output -- Prefix for output file.                   

        Output: 
        splicing_events_all_files_out -- List of filtered rMATS style 
            output file with splicing events from all input files. 
            File name- prefix + "_all_filtered_SE.MATS.JCEC.txt". 
        multiple_binding_splicing_out -- List of CSV files with splicing and
            multiple binding regions around a single spliced exon.
            The distance of the double binding regions was specified
            in the match_binding_to_splicing.py script. 
            File name- prefix + "_multiple_binding_splicing.csv".

    """
    # Uses argparse to organize command line arguments.
    args_to_parse = init_argparse()
    parsed = args_to_parse.parse_args()
    # Iterates through and executes main functions.
    find_output_all_splicing(parsed.filtered_splicing, parsed.output)
    find_output_multiple_binding(parsed.splicing_binding, parsed.output)

# Executes main function.
if __name__ == "__main__":
    main()