import argparse

def find_output_all_splicing(filtered_splicing, output):
    """Uses file list to find splicing events present in all
        files and outputs them in the same format.

        Arguments:
        filtered_splicing -- rMATS style output file list that is
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
        splicing_binding -- CSV file list with splicing and binding
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
                    #                   binding_region_start)
                    # NOTE- Since the splicing event does not contain the up-stream 
                    # region and the down-stream region, it could result in a false 
                    # positive. This is unlikely, but should be noted. This will likely 
                    # be changed in future revisions (revisions need to be performed 
                    # with match_splicing_to_binding.py). 
                    splicing_event = ("%s:%s:%s:%s:%s" % (line_list[2], line_list[3], line_list[4], 
                                                          line_list[5], line_list[9]))
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

def make_binding_list(splicing_binding):
    """Creates a dictionary for each file with the splicing event
        as the key and the binding region as the value.

        Arguments:
        splicing_binding -- CSV file list with splicing and binding
            information for each region. Produced by the 
            match_binding_to_splicing.py script.

        Output:
        binding_dictionary_list -- List of dictionaries with splicing 
            event as key and binding region as value. One is made 
            for each input file.
            {"Splicing Event":"clean_file_line"}.
        file_name_list -- List of file names. Order is used to coordinate
            with the binding_dictionary_list.
        title_line -- Title line for output.                  
    """
    # Iterates through each file and adds each splicing event
    # to a dictionary and then appands to a list.
    binding_dictionary_list = []
    file_name_list = []
    # Used to save the title line for output.
    title_line = ""
    for file in splicing_binding:
        binding_dictionary = {}
        # Collects file name.
        if file.find("/") != -1:
            file_name = file.split("/")[-1]
        # Adds file name to list.
        file_name_list.append(file_name)    
        with open(file, "r") as open_file:
            for line in open_file:
                line_clean = line.strip("\n")
                if line_clean[0] != "R":
                    line_list = line.split(",")
                    # splicing_event = (chromosome, exon-start, exon-end, strand,
                    #                   binding_region_start)
                    # NOTE- Since the splicing event does not contain the up-stream 
                    # region and the down-stream region, it could result in a false 
                    # positive. This is unlikely, but should be noted. This will likely 
                    # be changed in future revisions (revisions need to be performed 
                    # with match_splicing_to_binding.py). 
                    splicing_event = ("%s:%s:%s:%s:%s" % (line_list[2], line_list[3], line_list[4], 
                                                          line_list[5], line_list[9]))
                    binding_dictionary[splicing_event] = line_clean
                else:
                    title_line = (line_clean + ",Filename")  
        binding_dictionary_list.append(binding_dictionary)
    return binding_dictionary_list, file_name_list, title_line

def make_splicing_list(rmat_file):
    """Makes a list of dictionaries with all splicing events for
        each input file. Used to grab insignifcant splicing PSI 
        values if they are not present in the splicing_binding files.

        Arguments:
        rmat_file -- rMATS alternative skipped exon splicing 
            output file list (SE.MATS.JCEC.txt file). Should have 
            been generated using the same target organism 
            GTF file. All other types of alternative splicing 
            can be input except for retained intron (different 
            number of columns).

        Output:
        splicing_list -- List of dictionaries with all splicing events
            listed with filename. Used to grab insignifcant splicing
            PSI values as they would not be present in the 
            splicing_binding files.
            [{"Splicing Event:filename": PSI (Float)}].        
    """
    # Output list.
    splicing_list = []
    for file in rmat_file:
        # Gets filename. If directory name is included, removes it.
        if file.find("/") != -1:
            file_name = file.split("/")[-1]
        with open(file, "r") as open_file:
            # Used to save events per file.
            splicing_dictionary = {}
            for line in open_file:
                line_clean = line.strip("\n")
                if line_clean[0] != "I":
                    line_list = line.split("\t")
                    # Saves splicing event and PSI value.
                    # splicing_event = (chromosome, exon-start, exon-end, strand) 
                    splicing_event = ("%s:%s:%s:%s" % (line_list[3], line_list[5], line_list[6], 
                                                       line_list[4]))
                    splicing_dictionary[splicing_event] = float(line_list[22])
        splicing_list.append(splicing_dictionary)            
    return splicing_list

def make_normalized_count_list(normalized_counts, sample_csv_file):
    """Collects average counts for all binding regions. Used for 
        adding binding information for regions that are not signficant.

        Arguments:
        normalized_counts -- DESeq2 produced normalized counts file. 
            This file should contain all the normalized counts for all 
            protein concentrations.
        sample_csv_file -- Text file containing the file name to sample names in
            the normalized counts file. The filename should not contain the directory
            structure. This is used to match the files to the samples in the normalized
            counts file. The file is in CSV format. No title required. Used
            to grab counts that are not found in splicing_binding files.    

        Output:
        normalized_count_list -- List of dictionaries with all normalized counts
            for each input file.
            [{"binding_region": PSI (Float)}].    
    """
    # Grabs sample list to identify which columns to use for the average.
    sample_names = get_sample_names(sample_csv_file)
    # Grabs column locations for each sample.
    sample_location = get_sample_location(normalized_counts, sample_names)
    # Iterates through the normalized counts file and grabs the average
    # counts for each binding region.
    normalized_count_list = []
    # Iterates through file name in sample_location, saves each value,
    # calculates the average, and saves it to a dictionary.
    for file_name in sample_location:
        # Used to save the average counts.
        normalized_count_dictionary = {}
        # Opens normalized counts file.
        with open(normalized_counts, "r") as open_normalized_counts:
            # Used to save the average count for each binding region.
            for line in open_normalized_counts:
                # Filters for the title line.
                if line[1] != '"':
                    line_clean = line.strip("\n")
                    line_list = line_clean.split(",")
                    # Used to save the sum of the counts.
                    count_sum = 0
                    # Used to save the number of counts.
                    count_number = 0
                    # Iterates through the sample locations and adds the counts.
                    for location in sample_location[file_name]:
                        count_sum = count_sum + float(line_list[location])
                        count_number = count_number + 1
                    # Calculates the average.
                    average_count = (count_sum / count_number)
                    # Saves binding region name.
                    binding_region = line_list[0]
                    # Saves the average count to the dictionary.
                    normalized_count_dictionary[binding_region] = average_count
        normalized_count_list.append(normalized_count_dictionary)
    return normalized_count_list      
    
def get_sample_names(sample_csv_file):
    """Takes input sample_csv_file and returns a list of the sample names.

        Arguments:
        sample_csv_file -- Text file containing the file name to sample names in
            the normalized counts file. The filename should not contain the directory
            structure. This is used to match the files to the samples in the normalized
            counts file. The file is in CSV format. No title required.
            Format- filename,deseq_sample_name,..
            
        Output:    
        sample_names -- Dictionary of the sample names.
        {filename: [deseq_sample_name]}.
    """
    # Opens and reads sample_csv_file. Saves data to dictionary.
    sample_dictionary = {}
    # Count used for lines in file. Used to check length of 
    # sample_dictionary.
    line_count = 0
    with open(sample_csv_file) as open_sample_csv_file:
        for line in open_sample_csv_file:
            line_clean = line.strip("\n")
            line_list = line_clean.split(",")
            # Saves first column and dicionary key and
            # the remaining columns to a list since the
            # sample number can be variable.
            sample_dictionary[line_list[0]] = line_list[1:]
            # Checks to make sure content was added to the list.
            if len(line_list) > 1:
               line_count = line_count + 1      
    if (len(sample_dictionary) == 0 or line_count != len(sample_dictionary)):
        print("ERROR: No sample names found in sample_csv_file or sample"
              " names are not formatted correctly. System will continue"
              " to run, but samples might be missing. Please check file.")        
    return sample_dictionary

def get_sample_location(normalized_counts, sample_names):
    """Takes input normalized_counts file and returns a dictionary with the
        sample names as the key and the column location's list as the value.

        Arguments:
        normalized_counts -- DESeq2 produced normalized counts file. 
            This file should contain all the normalized counts for all 
            protein concentrations.
        sample_names -- Dictionary of the sample names.
            {filename: [deseq_sample_name]}.

        Output:    
        sample_location -- Dictionary of the sample names.
            {filename: [column_location]}.
    """
    # Dictionary used to save the file and column location for each sample.
    # {filename: [column_location]}
    sample_location = {}
    # Iterates through sample_names and uses the sample names to identify
    # which columns from the normalized_counts file to use for the average.
    for file_name in sample_names:
        # Used to save the sample column locations.
        sample_location_list = []
        sample_list = sample_names.get(file_name)
        # Opens normalized counts file and reads the first line.
        open_normalized_counts = open(normalized_counts, "r")
        first_line = open_normalized_counts.readline()
        open_normalized_counts.close()
        # Cleans and parses first line.
        first_line_list = first_line.strip("\n").split(",")
        # Iterates through each sample name and finds the corresponding
        # normalized counts. The column number is saved to a list.
        for sample in sample_list:
            # Iterates through each column in the first line.
            for column in first_line_list:
                # Checks to see if the sample names match. If so, saves the
                # column number.
                if column == sample:
                    sample_location_list.append(first_line_list.index(column))           
        sample_location[file_name] = sample_location_list
    return sample_location


def output_multiple_file_regions(binding_dictionary_list, splicing_list, title_line,
                                 normalized_count_list, sample_csv_file, 
                                 file_name_list, output):
    """Outputs splicing events that have been found in multiple files.

            Arguments:
            binding_dictionary_list -- List of dictionaries with splicing 
                event as key and binding region as value. One is made 
                for each input file.
                {"Splicing Event":"clean_file_line"}.
            splicing_list -- List of dictionaries with all splicing events
                listed with filename. Used to grab insignifcant splicing
                PSI values as they would not be present in the 
                splicing_binding files.
                {"Splicing Event:filename": PSI (Float)}.
            title_line -- Title line for output.        
            sample_csv_file -- Text file containing the file name to sample names in
                the normalized counts file. The filename should not contain the directory
                structure. This is used to match the files to the samples in the normalized
                counts file. The file is in CSV format. No title required. Used
                to grab counts that are not found in splicing_binding files
                Format- filename,deseq_sample_name,.
            file_name_list -- List of file names. Order is used to coordinate
                with the binding_dictionary_list.    
            normalized_counts -- DESeq2 produced normalized counts file. 
                This file should contain all the normalized counts for all 
                protein concentrations.
            output -- Prefix for output file. 

            Output:
            multiple_file_regions_out -- CSV file with splicing events
                and their respective binding regions that have been 
                found in multiple files. Multiple files are output,
                depending on the files they are found in. For example,
                The first file will output all splicing events in it
                that were also found in any subsequent file. The 
                distance of the double binding regions was specified
                in the match_binding_to_splicing.py script. 
                File name- prefix + "_multiple_file_regions.csv".
    """
    # Iterates through splicing/binding dictionaries, looks for splicing
    # events in later files, and outputs them to a file. Also pulls PSI
    # and average counts if the splicing event is not found in the
    # splicing_binding files.
    for index in range(len(binding_dictionary_list)):
        # Grabs file name from file_name_list based on
        # dictionary index number.
        file_name = file_name_list[index]
        # Opens output file for writing.
        multiple_file_regions_out = open(output + "_" + file_name + "_multiple_file_regions.csv", "w")
        # Writes title line to file.
        multiple_file_regions_out.write(title_line + "\n")
        # Iterates through splicing/binding dictionary.
        for splice_bind_event in binding_dictionary_list[index]:
            # Grabs splicing information only. Used to check if splicing
            # event is in other files. Only need first four variables.
            splicing_event_pre = splice_bind_event.split(":")
            splice_event_only = ":".join(splicing_event_pre[0:4])
            # List used to save common events until output.
            common_events = []
            # List used to save a index for each file. Used
            # to indicate for which files events have been found in
            # common.
            file_index_list = [index]
            # Used to save the splicing event and binding region.
            splicing_event_line = binding_dictionary_list[index][splice_bind_event]
            # Checks if splicing event is in any other files. Makes sure
            # to not check the same file. Also does not try to match lower
            # indices than the first index, as they have already been checked.
            for index2 in range(len(binding_dictionary_list)):
                if (index2 != index):
                    # Checks if the splicing event is in the dictionary.
                    if splice_bind_event in binding_dictionary_list[index2]:
                        # Grabs file name from file_name_list based on
                        # dictionary index number for file 2.
                        file_name_2 = file_name_list[index2]
                        # If splicing event is in dictionary, adds line to list.
                        splicing_event_line_next = binding_dictionary_list[index2][splice_bind_event]
                        # Saves binding region name if the counts from the normalized counts
                        # file is required to be collected.
                        binding_region = splicing_event_line_next.split(",")[0]
                        # Adds index to list.
                        file_index_list.append(index2)
                        # Adds splicing events to list. If first entry, adds from both files.
                        if len(common_events) == 0:
                            common_events.append(splicing_event_line + "," + file_name)
                            common_events.append(splicing_event_line_next + "," + file_name_2)
                        else:
                            common_events.append(splicing_event_line_next + "," + file_name_2) 
            # Counts the number of files the splicing event was found in.
            if len(common_events) > 1:
                # Checks for any missing splicing events in the splicing_binding files.
                for index_check in range(len(binding_dictionary_list)):
                    if index_check not in file_index_list:
                        # Grabs file name from file_name_list based on
                        # dictionary index number.
                        file_name_check = file_name_list[index_check]
                        # Checks if splicing event is in the respective dictionary insplicing_list.
                        if splice_event_only in splicing_list[index_check]:
                            # Grabs PSI value.
                            psi_value = splicing_list[index_check][splice_event_only]                  
                        # Gets average counts for the binding region.
                        if binding_region in normalized_count_list[index_check]:
                            average_count = normalized_count_list[index_check][binding_region]
                        # Saves splicing event and binding region to a string, with all
                        # other values as "NA".
                        missing_splicing_event = (binding_region + "," + "".join(["NA,"] * 6) 
                                                  + str(psi_value) + "," + "".join(["NA,"] * 8) 
                                                  + str(average_count) + "," + file_name_check)
                        # Adds missing splicing event to common events.
                        common_events.append(missing_splicing_event)
                # Outputs splicing events to file.
                for line in common_events:
                    multiple_file_regions_out.write("%s\n" % (line))
        # Closes output file.
        multiple_file_regions_out.close()
      
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
                        help="rMATS style output file list that is \
                              filtered by delta PSI and an FDR of 10 percent or less. \
                              Contains the 'filtered_SE.MATS.JCEC.txt' suffix.",                          
                        required = True)
    parser.add_argument("-b", "--splicing_binding", action = "store", type = str, nargs='+',
                        help="CSV file list with splicing and binding \
                              information for each region. Produced by the \
                              match_binding_to_splicing.py script.",
                        required = True)
    parser.add_argument("-c", "--sample_csv_file", action = "store", type = str,
                        help="Text file containing the file name to sample names in \
                              the normalized counts file. The filename should not contain the directory \
                              structure. This is used to match the files to the samples in the normalized \
                              counts file. The file is in CSV format. No title required. Used \
                              to grab counts that are not found in splicing_binding files \
                              Format- filename,deseq_sample_name,.",
                        required = True)
    parser.add_argument("-r", "--rmats_file", action = "store", type = str, nargs='+',
                        help="rMATS alternative skipped exon splicing \
                              output file list (SE.MATS.JCEC.txt file). Should have \
                              been generated using the same target organism \
                              GTF file. All other types of alternative splicing \
                              can be input except for retained intron (different \
                              number of columns).",
                        required = True)
    parser.add_argument("-n", "--normalized_counts", action = "store", type = str,
                        help="DESeq2 produced normalized counts file. \
                              This file should contain all the normalized counts for all \
                              protein concentrations.",
                        required = True)
    parser.add_argument("-o", "--output", action = "store", type = str,
                        help="Prefix for output files.",
                        required = True)                                                                 
    return parser

def main():
    """Main script function. Executes other functions.

        Arguments:
        filtered_splicing -- rMATS style output file list that is
            filtered by delta PSI and an FDR of 10% or less.
            Contains the 'filtered_SE.MATS.JCEC.txt' suffix.
        splicing_binding -- CSV file list with splicing and binding
            information for each region. Produced by the 
            match_binding_to_splicing.py script. Order of files in 
            the file should be the same as the order of the files
            in other input files.
        sample_csv_file -- Text file containing the file name to sample names in
            the normalized counts file. The filename should not contain the directory
            structure. This is used to match the files to the samples in the normalized
            counts file. The file is in CSV format. No title required. Used
            to grab counts that are not found in splicing_binding files. Order
            of files in the file should be the same as the order of the files
            in other input files.
            Format- filename,deseq_sample_name,.
        rmat_file -- rMATS alternative skipped exon splicing 
            output file list (SE.MATS.JCEC.txt file). Should have 
            been generated using the same target organism 
            GTF file. All other types of alternative splicing 
            can be input except for retained intron (different 
            number of columns). Order of files in the file should 
            be the same as the order of the files in other input files.   
        normalized_counts -- DESeq2 produced normalized counts file. 
            This file should contain all the normalized counts for all 
            protein concentrations. Used to grab counts that are not found in
            splicing_binding files.           
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
        multiple_file_regions_out -- CSV file with splicing events
            and their respective binding regions that have been 
            found in multiple files. Multiple files are output,
            depending on the files they are found in. For example,
            The first file will output all splicing events in it
            that were also found in any subsequent file. The 
            distance of the double binding regions was specified
            in the match_binding_to_splicing.py script. 
            File name- prefix + "_multiple_file_regions.csv".    
    """
    # Uses argparse to organize command line arguments.
    args_to_parse = init_argparse()
    parsed = args_to_parse.parse_args()
    # Iterates through and executes main functions.
    find_output_all_splicing(parsed.filtered_splicing, parsed.output)
    find_output_multiple_binding(parsed.splicing_binding, parsed.output)
    binding_dictionary_list, file_name_list, title_line = make_binding_list(parsed.splicing_binding)
    splicing_list = make_splicing_list(parsed.rmats_file)
    normalized_count_list = make_normalized_count_list(parsed.normalized_counts, parsed.sample_csv_file)
    output_multiple_file_regions(binding_dictionary_list, splicing_list, title_line,
                                 normalized_count_list, parsed.sample_csv_file, 
                                 file_name_list, parsed.output)

# Executes main function.
if __name__ == "__main__":
    main()