import argparse

def process_regions(deseq_file_list, normalized_counts, sample_csv_file, output):
    """Takes input DeSeq2 produced files (normalized counts or differential
        expression) and executes the data through other functions to output
        the normalized counts file with the lower concentration binding regions
        removed.

        Arguments:
        deseq_file_list -- DESeq2 produced file type. Should contain all files for 
            all protein concentrations. All the files shoule be the SM filtered
            and filtered with the control sample (i.e., 
            x50_vs_uni_significant_filt_annot.csv). The files should be listed 
            in the order desired to be filtered. This is is usually the lower 
            concentration to higher concentration. The file is in CSV format.
        normalized_counts -- DESeq2 produced normalized counts file. This file
            should contain all the normalized counts for all protein concentrations.
        sample_csv_file -- Text file containing the file name to sample names in
            the normalized counts file. The filename should not contain the directory
            structure. This is used to match the files to the samples in the normalized
            counts file. The file is in CSV format. No title required.
            Format- filename,deseq_sample_name,..
        output -- Directory output for output file.
               
        Output:    
        Uses output_file function to output the normalized counts file with the
            lower concentration binding regions removed. The file is in CSV format.
    """
    # Gets the sample names from the sample_csv_file.
    sample_names = get_sample_names(sample_csv_file)
    # Uses the sample names to get the average normalized counts 
    # for each sample.
    average_counts = get_average_counts(normalized_counts, sample_names)
    # Used to remove the lower concentration binding regions.
    lower_region, all_regions, title_dictionary  = remove_lower_regions(deseq_file_list)
    # Outputs the region files with the average_counts for each region
    # added to the end of each line.
    output_file(average_counts, lower_region, all_regions, title_dictionary, output)
   

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

def get_average_counts(normalized_counts, sample_names):
    """Takes input normalized_counts file and returns a list of the average
        normalized counts for each sample file.

        Arguments:
        normalized_counts -- DESeq2 produced normalized counts file. 
            This file should contain all the normalized counts for all 
            protein concentrations.
        sample_names -- Dictionary of the sample names.
            {filename: [deseq_sample_name]}.

        Output:    
        average_counts -- Dictionary of the average normalized counts 
            for each sample.
            {sample_name: {binding_region: average_normalized_count}}.
    """
    # Grabs sample column location in normalized_counts file. Returns
    # a dictionary with the following format.
    # {filename: [column_location]}
    sample_location = get_sample_location(normalized_counts, sample_names)
    # Dictionary used to save the average normalized counts for each sample.
    # {file: {binding_region: average_normalized_count}}
    average_counts = {}
    # Iterates through sample_location and uses the sample names to identify
    # which columns from the normalized_counts file to use for the average.
    for file_name in sample_location:
        average_counts[file_name] = {}
        # Used to save the sample column locations.
        sample_location_list = sample_location.get(file_name)
        # Opens normalized counts file and loops through each line.
        # Uses the sample_location_list to grab the columns specific
        # to the different samples for a single protein concentration.
        with open(normalized_counts) as open_normalized_counts:
            # Dictionary used to save counts until appended to average_counts.
            # {binding_region: [normalized_counts]}
            region_normalized_counts = {}
            for line in open_normalized_counts:
                line_list = line.strip("\n").split(",")
                # Uses the sample_location_list to grab the columns specific
                # to the different samples for a single protein concentration.
                # If statement skips title line.
                if len(line_list[0]) > 1:
                    binding_region = line_list[0]
                    # Saves the normalized counts for each sample.
                    normalized_counts_list = []
                    for location in sample_location_list:
                        normalized_counts_list.append(line_list[location])     
                    # Grabs lenth of normalized_counts_list to use for average.
                    length = len(normalized_counts_list)
                    # Count used for summing values.
                    count_sum = 0
                    for count in normalized_counts_list:
                        count_sum = (count_sum + float(count))
                    # Calculates average. If length is 0, average is 0.       
                    average = (count_sum / length) if length > 0 else 0
                    # Grabs dictionary for the file and adds the binding region
                    # and average normalized counts.
                    region_normalized_counts = average_counts.get(file_name)    
                    # Saves the binding region and the normalized counts for region.
                    region_normalized_counts[binding_region] = average
                    # Saves the sample name and the normalized counts for the sample file.
                    average_counts[file_name] = region_normalized_counts             
    return average_counts
          

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

def remove_lower_regions(deseq_file_list):
    """Takes in a list of DESeq2 filtered files and removes any binding
        region that was found in a file previously in the list. The
        list should be in the order of lower to higher concentration.

        Arguments:
        deseq_file_list -- DESeq2 produced file type. Should contain all files for 
            all protein concentrations. All the files shoule be the SM filtered
            and filtered with the control sample (i.e., 
            x50_vs_uni_significant_filt_annot.csv). The files should be listed 
            in the order desired to be filtered. This is is usually the lower 
            concentration to higher concentration. The file is in CSV format.

        Output:
        lower_regions -- Dictionary of binding regions that were found in the
            lower concentration files.
            {"file_name": {binding_region: line}}. 
        all_regions -- Dictionary of all binding regions found in the files.
            {"file_name": {binding_region: line}}.
        title_dictionary -- Dictionary of the title lines for each file.
            Used to output the title lines to the output files.
            {"file_name": title_line}.    
    """
    # Dictionaries used to store binding regions per sample file and the title
    # lines. These dictionaries are returned from the function.
    lower_region = {}
    all_regions = {}
    title_dictionary = {}
    # Set used to store all lower regions.
    lower_region_set = set()
    # Iterates through each file in the list.
    for file in deseq_file_list:
        # Saves the file name to a string.
        file_name = ""
        if file.find("/") != -1:
            file_name = file.split("/")[-1]
        else:
            file_name = file
        # Adds file name to dictionaries for easier saving later.    
        lower_region[file_name] = {}
        all_regions[file_name] = {}   
        # Opens file and reads each line.
        with open(file) as open_file:
            for line in open_file:
                line_clean = line.strip("\n")
                line_list = line_clean.split(",")
                # Checks for title line. Depends on length of 
                # string in the first column.
                if len(line_list[0]) > 2:
                    # Checks to see if the binding region has already been found.
                    # If so, saves the binding region to the lower_region dictionary.
                    # If not, saves the binding region to the lower_region_set.
                    if line_list[0] not in lower_region_set:
                        #lower_region_temp[line_list[0]] = line_clean
                        lower_region[file_name][line_list[0]] = line_clean
                        lower_region_set.add(line_list[0])   
                    all_regions[file_name][line_list[0]] = line_clean
                else:
                    title_dictionary[file_name] = line_clean
    # Printing only keys from dictionaries for testing.                               
    return lower_region, all_regions, title_dictionary             

def output_file(average_counts, lower_region, all_regions, 
                title_dictionary, output):
    """Takes dictionaries with the average normalized counts,
        one with the lower regions filtered out, and one with
        all regions and outputs the data into two files.

        Arguments:
        average_counts -- Dictionary of the average normalized counts 
            for each sample.
            {sample_name: {binding_region: average_normalized_count}}
        lower_region -- Dictionary of binding regions that were found in the
            lower concentration files.
            {"file": {binding_region: line}} .
        all_regions -- Dictionary of all binding regions found in the files.
            {"file": {binding_region: line}}.
        title_dictionary -- Dictionary of the title lines for each file.
            Used to output the title lines to the output files.
            {"file": title_line}.
        output -- Directory output for output file.    

        Output:
        deseq_file_out -- Same format as deseq_file, but with one extra column named
            "normalized_counts". This column contains the average DESeq2 normalized
            counts for the sample set. The file contains only the binding regions
            that were first bound at this concentration. This means the lower concentration
            significant binding regions have been removed. The file is in CSV format and
            will contain the file name with "no_lower" added to the end of the filename
            ".csv" suffix.
        normalized_counts_out -- Same as the input DESeq2 file but with an extra column
            added to the end. This column contains the average DESeq2 normalized counts
            for the sample set. The file is in CSV format and will contain the file name
            with "counts" added to the end of the filename ".csv" suffix.        
    """
    # Iterates through all_regions and outputs the data to two files. Only one loop is needed
    # since all_regions contains all the binding regions for each file.
    for file_name in all_regions:
        # Removes the ".csv" suffix from the file name.
        file_name_out = file_name.replace(".csv", "")
        # Opens output files and writes the title line. Only used for filename
        # output. The title line is not used in the output files.
        if output[-1] == "/":
            deseq_file_out = open(output + file_name_out + "_no_lower.csv", "w")
            normalized_counts_out = open(output + file_name_out + "_counts.csv", "w")
        else:
            deseq_file_out = open(output + "/" + file_name_out + "_no_lower.csv", "w")
            normalized_counts_out = open(output + "/" + file_name_out + "_counts.csv", "w")       
        deseq_file_out.write(title_dictionary.get(file_name) + ",normalized_counts\n")
        normalized_counts_out.write(title_dictionary.get(file_name) + ",normalized_counts\n")    
        # Iterates through each binding region in the file.
        for region in all_regions.get(file_name):
            # Grabs normalized counts for the region.
            normalized_counts = average_counts.get(file_name).get(region)
            # Writes the binding region and normalized counts to the output files.
            if region in lower_region.get(file_name):
                deseq_file_out.write(lower_region.get(file_name).get(region) + "," + str(normalized_counts) + "\n")
            normalized_counts_out.write(all_regions.get(file_name).get(region) + "," + str(normalized_counts) + "\n")
        deseq_file_out.close()
        normalized_counts_out.close()
        
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
    parser.add_argument("-de", "--deseq_file_list", action = "store", type = str, nargs='+',
                        help="DESeq2 produced file type. Should contain all files for \
                              all protein concentrations. All the files shoule be the SM filtered \
                              and filtered with the control sample (i.e., \
                              x50_vs_uni_significant_filt_annot.csv). The files should be listed \
                              in the order desired to be filtered. This is is usually the lower \
                              concentration to higher concentration. The file is in CSV format.", 
                              required = True)
    parser.add_argument("-n", "--normalized_counts", action = "store", type = str, 
                        help = "DESeq2 produced normalized counts file. This file \
                                should contain all the normalized counts for all \
                                protein concentrations..", required = True)
    parser.add_argument("-csv", "--sample_csv_file", action = "store", type = str,
                        help = "Text file containing the file name to sample names in \
                                the normalized counts file. The filename should not contain \
                                the directory structure. This is used to match the files to \
                                the samples in the normalized counts file. The file is in CSV \
                                format. No title required. \
                                Format- filename,deseq_sample_name,...", required = True)
    parser.add_argument("-o", "--output", action = "store", type = str, default = "./",
                        help = "Directory output to add to the output file name.", 
                                required = False)                                                                          
    return parser

def main():
    """Main script function. Executes other functions.

        Arguments:
        deseq_file_list -- DESeq2 produced file type. Should contain all files for 
            all protein concentrations. All the files shoule be the SM filtered
            and filtered with the control sample (i.e., 
            x50_vs_uni_significant_filt_annot.csv). The files should be listed 
            in the order desired to be filtered. This is is usually the lower 
            concentration to higher concentration. The file is in CSV format.
        normalized_counts -- DESeq2 produced normalized counts file. This file
            should contain all the normalized counts for all protein concentrations.
        sample_csv_file -- Text file containing the file name to sample names in
            the normalized counts file. The filename should not contain the directory
            structure. This is used to match the files to the samples in the normalized
            counts file. The file is in CSV format. No title required.
            Format- filename,deseq_sample_name,..
        output -- Directory output for output file.   
            
        Output:    
        deseq_file_out -- Same format as deseq_file, but with one extra column named
            "normalized_counts". This column contains the average DESeq2 normalized
            counts for the sample set. The file contains only the binding regions
            that were first bound at this concentration. This means the lower concentration
            significant binding regions have been removed. The file is in CSV format and
            will contain the file name with "no_lower" added to the end of the filename
            ".csv" suffix.
        normalized_counts_out -- Same as the input DESeq2 file but with an extra column
            added to the end. This column contains the average DESeq2 normalized counts
            for the sample set. The file is in CSV format and will contain the file name
            with "counts" added to the end of the filename ".csv" suffix.     
    """
    # Uses argparse to organize command line arguments.
    args_to_parse = init_argparse()
    parsed = args_to_parse.parse_args()
    # Iterates through and executes main functions.
    process_regions(parsed.deseq_file_list, parsed.normalized_counts,
                    parsed.sample_csv_file, parsed.output)

# Executes main function.
if __name__ == "__main__":
    main()