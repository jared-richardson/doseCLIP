import argparse

def process_files(deseq_file_list, output_prefix, output):
    """Iterates through each file, adds normalized counts to a dictionary,
        and then outputs files.

    Arguments:
    deseq_file_list -- DESeq2 produced file type. Should contain all files for 
        all protein concentrations. All the files shoule be the SM filtered
        and filtered with the control sample (i.e., 
        x50_vs_uni_significant_filt_annot.csv). The files should be listed 
        in the order desired to be filtered. This is is usually the lower 
        concentration to higher concentration. The file is in CSV format.
        Must contain last average counts from remove_lower_regions.py script.
        Should contain an entire file set for each protein concentration. Files
        with lower concentration regions should be executed separately from
        files without them.
    output_prefix -- Prefix for output files.
    output -- Directory output for output file.

    Output:
    normalized_average_counts_out -- File contains average normalized counts for each
        sample file based on file input. File is in CSV format.
        File name- "output" + "sample_name" + "_average_normalized_counts.csv"
        Format- "Sample_name,Average_normalized_counts\n" 
    all_normalized_counts_out -- File contains all normalized counts for each
        single sample file based on file input. File is in CSV format.
        File name- "output" + "sample_name" + "_average_normalized_counts.csv"
        Format- Format is file name with counts below.
    all_normalized_counts_out_combined -- File contains all normalized counts for each
        single sample file's counts from another file. File is in CSV format.
        File name- "output" + "sample_name" + "_average_normalized_counts.csv"
        Format- Format is file name with counts below.     
    """
    # Opens average output file and writes the title line. Only used for filename
    # output. The title line is not used in the output files.
    if output[-1] == "/":
        normalized_average_counts_out = open(output + output_prefix + "_average_normalized_counts.csv", "w")
    else:
        normalized_average_counts_out = open(output + "/" + output_prefix + "_average_normalized_counts.csv", "w")
    # Writes out title for normalized_average_counts_out.
    normalized_average_counts_out.write("Sample_name,Average_normalized_counts\n")    
    # Intitializes dictionary to store each sample file's normalized counts.
    sample_dictionary = {}
    # Average count dictionary. Stores the average normalized counts for each sample.
    average_count_dictionary = {}
    # Iterates through each DESeq2 file and adds the normalized counts to a dictionary.
    for file in deseq_file_list:
        count_dictionary, sample_name = get_normalized_counts(file)
        sample_dictionary[sample_name] = count_dictionary    
    # Gets average count for each sample. Also finds binding regions that are present
    # in each other sample's dictionary and gets average if they are present. Saves
    # the name of the combination of samples as sample_name-sample_name.
    # Also outputs each count to a CSV file in case individual counts are needed.        
    for sample in sample_dictionary:
        # Opens sample output file for writing.
        if output[-1] == "/":
            all_normalized_counts_out = open(output + sample + "_all_normalized_counts.csv", "w")
        else:
            all_normalized_counts_out = open(output + "/" + sample + "_all_normalized_counts.csv", "w")
        all_normalized_counts_out.write(sample + "\n")    
        average_count = 0
        count = 0
        # Count used so that sample average is only calculated once.
        sample_flag_count = 0
        for other_sample in sample_dictionary:
            if sample != other_sample:
                average_count_both = 0
                count_both = 0
                sample_flag_count += 1
                # Combined sample name for output.
                combined_name = (sample + "-" + other_sample)
                # Opens sample output file for writing.
                if output[-1] == "/":
                    all_normalized_counts_out_combined = open(output + combined_name + "_all_normalized_counts.csv", "w")
                else:
                    all_normalized_counts_out_combined = open(output + "/" + combined_name + "_all_normalized_counts.csv", "w")
                all_normalized_counts_out_combined.write(combined_name + "\n")
                for region in sample_dictionary.get(sample):
                    if sample_flag_count == 1:
                        # Filters for outliers.
                        if sample_dictionary.get(sample).get(region) < 4000:
                            average_count += sample_dictionary.get(sample).get(region)
                            count += 1
                            # Writes out all normalized counts to a file.
                            all_normalized_counts_out.write(str(sample_dictionary.get(sample).get(region)) + "\n")
                    if region in sample_dictionary.get(other_sample):
                        # Filters for outliers.
                        if sample_dictionary.get(sample).get(region) < 4000:
                            average_count_both += sample_dictionary.get(other_sample).get(region)
                            count_both += 1
                            # Writes out all normalized counts to a file.
                            all_normalized_counts_out_combined.write(str(sample_dictionary.get(other_sample).get(region)) + "\n")
                if count_both != 0:        
                    average_count_dictionary[combined_name] = (average_count_both / count_both)             
        average_count_dictionary[sample] = (average_count / count)
        all_normalized_counts_out.close()
        all_normalized_counts_out_combined.close()
    # Writes out average counts to a file.
    for sample in average_count_dictionary:
            normalized_average_counts_out.write(sample + "," + str(average_count_dictionary.get(sample)) + "\n")
    normalized_average_counts_out.close()    

def get_normalized_counts(deseq_file):
    """Adds normalized counts to a dictionary.

    Arguments:
    deseq_file -- DESeq2 produced file type. Should contain all files for 
        all protein concentrations. All the files shoule be the SM filtered
        and filtered with the control sample (i.e., 
        x50_vs_uni_significant_filt_annot.csv). The file is in CSV format.
        Must contain last average counts from remove_lower_regions.py script.
        Should contain an entire file set for each protein concentration. Files
        with lower concentration regions should be executed separately from
        files without them.

    Output:
    count_dictionary -- Dictionary of the average normalized counts 
        for each sample.
        {binding_region: normalized_count}
    sample_name -- Name of the sample file. Used for output file name.       
    """
    # Intializes output dictionary.
    count_dictionary = {}
    # Gets filename for to save as sample name.
    if deseq_file.find("/") != -1:
        sample_name = deseq_file.split("/")[-1].strip(".csv")
    else:
        sample_name = deseq_file.strip(".csv")            
    # Opens file and reads each line.
    with open(deseq_file) as open_file:
        for line in open_file:
            line_clean = line.strip("\n")
            line_list = line_clean.split(",")
            # Checks for title line. Depends on length of 
            # string in the first column.
            if len(line_list[0]) > 2:
                # Adds normalized counts to dictionary.
                count_dictionary [line_list[0]] = float(line_list[-1])
    return (count_dictionary, sample_name)
        
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
                              concentration to higher concentration. The file is in CSV format. \
                              Must contain last average counts from remove_lower_regions.py script. \
                              Should contain an entire file set for each protein concentration. Files \
                              with lower concentration regions should be executed separately from \
                              files without them.", 
                              required = True)
    parser.add_argument("-p", "--output_prefix", action = "store", type = str,
                        help="Prefix for output files.",
                        required = True)
    parser.add_argument("-o", "--output", action = "store", type = str, default = "./",
                        help = "Directory to output files.", 
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
            Must contain last average counts from remove_lower_regions.py script.
            Should contain an entire file set for each protein concentration. Files
            with lower concentration regions should be executed separately from
            files without them.
        output_prefix -- Prefix for output files.
        output -- Directory output for output file.   
            
        Output:    
        normalized_average_counts_out -- File contains average normalized counts for each
        sample file based on file input. File is in CSV format. 
        all_normalized_counts_out -- File contains all normalized counts for each
        sample file based on file input. File is in CSV format.   
    """
    # Uses argparse to organize command line arguments.
    args_to_parse = init_argparse()
    parsed = args_to_parse.parse_args()
    # Iterates through and executes main functions.
    process_files(parsed.deseq_file_list, parsed.output_prefix, parsed.output)

# Executes main function.
if __name__ == "__main__":
    main()