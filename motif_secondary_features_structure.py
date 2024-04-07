import argparse

def organize_nucleotides(fasta_file_name, kmer_size, motif_list = ["YGCY"], 
                         motif_title = "YGCY"):
    """Analyzes a region FASTA file for the presence of a motif.

        Arguments:
        fasta_file_name -- String, name of FASTA file being analyzed.
        kmer_size -- Size of kmer to be analyzed.
        motif_list -- List of motifs to be analyzed. Default is YGCY.
        motif_title -- Title for motif output file. Default is YGCY.

        Output:    
        region_dictionary -- Dictionary of nucleotide lists for the nucleotides
            surrounding each motif and the motifs themselves. The data is arranged
            in a list of four lists. The descriptions of each list are listed below.
            "nucleotide_list"- List, string list of nucleotides, seven nucleotides
                around all motifs.
            "nucleotide_intra_list"- List, string list of nucleotides between motifs
                in a motif group. A motif group is defined as motifs within seven 
                nucleotides of each other.
            "nucleotide_inter_list"- List, string list of nucleotides between motifs
                separated by more than seven nucleotides.
            "motifs_group_list"- List, string list of motifs and surrounding 
                nucleotides within a group.
            region_size -- Integer, size of binding region (last variable).    
            {"region_name": [nucleotide_list, nucleotide_intra_list,
                             nucleotide_inter_list, motifs_group_list,
                             region_size].
        structure_dictionary -- Dictionary of secondary structure predictions for
            motifs and surrounding nucleotides. The data is a dictionary with the
            region as the key and the list of secondary structure predictions as
            the value, arranged in three lists.
            {"region_name": [structure, y_structure, gc_structure]}.                       

    """
    # Used to store nucleotide data. See description above.
    region_dictionary = {}
    structure_dictionary = {}
    # Checks motif_list for "YGCY", the default. If found, replaces
    # motif_list populates with all YGCY variants. If not found,
    # the input motif list is used. Uses motif_title to name output.
    if motif_list == ["YGCY"]:
         motif_list = ["CGCC", "CGCU", "UGCC", "UGCU"]
         motif_title = "YGCY"
    else:
         motif_list = motif_list
         motif_title = motif_title
    # Used to store name of binding region.
    region_name = ""
    open_file = open(fasta_file_name, 'r')
    # Use to detect triplets of lines that refer to the
    # same region. 0 is reset, and 1 is stay the same.
    line_count = 0
    # Loops through each line of the FASTA file. If title line
    # title is saved, if nucleotides, motif counts are generated.        
    for line in open_file:
        line_clean = line.strip("\n")
        # Only resets data structures if the line is a title line.
        if line_count == 0:
            line_count = 1
            # Used to grab region size.
            region_size = 0
            # List, string list of nucleotides, seven nucleotides 
            # around all motifs. String is used to store nucleotides
            # for each group.
            nucleotide_list = []
            nucleotide_string = ""
            # List, string list of nucleotides between motifs
            # in a motif group. A motif group is defined as 
            # motifs within seven nucleotides of each other.
            # String is used to store nucleotides for each group.
            nucleotide_intra_list = []
            nucleotide_intra_string = ""
            # Saves index of nucleotide intra list.
            nucleotide_intra_index = []
            # List, string list of nucleotides between motifs
            # separated by more than seven nucleotides.
            # String is used to store nucleotides for each group.
            nucleotide_inter_list = []
            nucleotide_inter_string = ""
            # List, string list of motifs and surrounding 
            # nucleotides within a group. String is used to
            # store nucleotides for each group. Groups need
            # to be separted with this list so the groups can
            # be counted and analyzed separetly.
            motifs_group_list = []
            motifs_group_string = ""
            # Lists used to store indexes of nucleotides around
            # motifs (structure_index), indexes of "Y" in "YGCY"
            # (y_index), and indexes of "GC" in "YGCY" (gc_index).
            structure_index = []
            y_index = []
            gc_index = []
            # Lists used to store secondary structure predictions
            # for each nucleotide set. Same order as above.
            structure = []
            y_structure = []
            gc_structure = []
            # Used to detect nucleotide line. Every third line is
            # structure, so only lines with nucleotides are analyzed.
        if (line[0] != ">") and (line[0] in ["A", "C", "G", "U"]):
            # Grabs region size.
            region_size = len(line_clean)
            # Nucleotide and structure lists. List used to keep 
            # track of nucleotides and structure around
            # motif of interest. Does not contain motifs.
            around_motif = []
            around_index = []
            # Distance between each motif.
            motif_spacing = 0
            # Initializes motif_count for each kmer.
            motif_count = 0
            # Used to determine if a motif has been found previously in last
            # four nucleotides. If so changes the spacing and motif counts.
            double_motif_count = 0
            # Count used for the last part of the spacing between motifs.
            # If a motif occurs that has overlap in the seven nucleotide
            # spacer, the counter will monitor this and determine the
            # next nucleotides were part of the same motif group.
            motif_counter = 0
            # Count used to grab the last motif group data if the space between
            # motifs is large.
            large_space_count = 0
            # Flag used to determine if the last motif was a double motif. If so,
            # adds the last three nucleotides to the nucleotide string, and only
            # three nucleotides of the current kmer to the motifs_group_string.
            double_motif_count_flag = False
            # Loops through each nucleotide line, identifies wanted kmers,
            # and organizes the nucleotides into lists.
            for start in range(len(line_clean) - kmer_size + 1):     
                kmer = line_clean[start:start + kmer_size]
                # Adds index of nucleotides around motif to around_index list.
                kmer_indexes = [start + index for index in range(kmer_size)]
                # Adds nucleotides to around_motif list. If list is empty,
                # adds all nucleotides. If not, adds only the last nucleotide.
                if (motif_count == 0 and len(around_motif) == 0):
                    for nucleotide in kmer:
                        around_motif.append(nucleotide)
                    for index in kmer_indexes:
                        around_index.append(index)
                    motif_spacing += 4
                else:
                    around_motif.append(kmer[-1])
                    around_index.append(kmer_indexes[-1])
                    motif_spacing += 1
                # If kmer is in motif_list, saves up to the last seven nucleotides
                # and appends to nucleotide_list. Also determines what other
                # nucleotides need to be appended to other lists.
                if kmer in motif_list:
                    # Adds index of motifs to their respecctive lists.
                    # The first character and the last character of the motif
                    # index is added to y_index and the second and third
                    # characters are added to gc_index.
                    y_index.append(kmer_indexes[0])
                    y_index.append(kmer_indexes[-1])
                    gc_index.append(kmer_indexes[1])
                    gc_index.append(kmer_indexes[2])
                    # Removes last three nucleotides from around_motif list if
                    # there was a double motif (i.e., YGCYGCY). If not, removes
                    # the standard four nucleotides. Set double_motif_count_flag
                    # to true so only three nucleotides are added to motif_group_string.
                    if double_motif_count != 0:
                        del around_motif[-3:]
                        del around_index[-3:]   
                        motif_spacing -= 3
                        double_motif_count_flag = True 
                    else:
                        del around_motif[-4:]
                        del around_index[-4:]
                        motif_spacing -= 4    
                    # Used for first motif in binding region.
                    if motif_count == 0:
                        # Grabs previous seven non-motif nucleotides
                        # and adds them to nucleotide strings.
                        last_nucleotides = around_motif[-7:]
                        last_index = around_index[-7:]
                        for nucleotide in last_nucleotides:
                            nucleotide_string += nucleotide
                            motifs_group_string += nucleotide
                        for index in last_index:
                            structure_index.append(index)   
                    # If previous motif has been found, checks for the 
                    # double_motif_count_flag. If true no nucleotides
                    # are added since they would have already been added.        
                    elif double_motif_count_flag != True:
                        # Checks to see if motif is within seven nucleotides,
                        # this means that the motif is in the same group. If
                        # not, adds nucleotides to nucleotide_inter_list.
                        last_nucleotides = around_motif[-motif_spacing:]
                        last_index = around_index[-motif_spacing:]  
                        if (0 < motif_spacing <= 7):
                            # Checks to see if motif was actually within seven
                            # nucleotides using large_space_count. If so, 
                            # grabs last item from motifs_group_list, deletes
                            # this item, and removes the last seven nucleotides
                            # from motifs_group_string and nucleotide_string, this
                            # way the nucleotides are not added twice.
                            if large_space_count > 0:
                                motifs_group_string = motifs_group_list[-1]
                                del motifs_group_list[-1]
                                motifs_group_string = motifs_group_string[:-7]
                                nucleotide_string = nucleotide_string[:-7]
                                structure_index = structure_index[:-7]
                            # Adds to nucleotides to applicable strings.
                            for nucleotide in last_nucleotides:
                                nucleotide_string += nucleotide
                                motifs_group_string += nucleotide
                                nucleotide_intra_string += nucleotide    
                            for index in last_index:
                                structure_index.append(index)
                                nucleotide_intra_index.append(index) 
                        # Makes sure there are actually nucleotides to add and
                        # then adds the nucleotides to the appropriate lists.          
                        elif (motif_spacing != 0):
                            # Adds to all nucleotides to nucleotide_inter_string.
                            for nucleotide in last_nucleotides:
                                nucleotide_inter_string += nucleotide
                            # Checks to make sure there the distance was actually
                            # greater than seven. If not, removes last added
                            # nucleotides and adds current nucleotides to the
                            # required lists. nucleotide_intra_string is
                            # added to nucleotide_string because this ensures
                            # no double additions after the deletion.
                            if large_space_count > 0:
                                motifs_group_string = motifs_group_list[-1]
                                del motifs_group_list[-1]
                                nucleotide_string = nucleotide_string[:-7]
                                nucleotide_string += nucleotide_intra_string
                                around_index = around_index[:-7]
                                around_index += nucleotide_intra_index
                                # Resets count bec  ause a motif was found.
                                large_space_count = 0   
                            # Now adds the last seven nucleotieds to the
                            # appropriate lists.    
                            last_nucleotides = around_motif[-7:]
                            last_index = around_index[-7:]
                            for nucleotide in last_nucleotides:
                                nucleotide_string += nucleotide
                                motifs_group_string += nucleotide
                            for index in last_index:
                                structure_index.append(index)                                                  
                    motif_spacing = 0
                    motif_count += 1
                    # Ensures the right motif nucleotides are a
                    # added if a double motif was found.
                    if double_motif_count_flag == True:
                        kmer = kmer[-3:]
                        motifs_group_string += kmer
                    else:    
                        motifs_group_string += kmer
                    # Resets double_motif_count_flag and sets double_motif_count
                    # to three to check for future double motifs.
                    double_motif_count_flag = False
                    double_motif_count = 3
                # Performs the correct procedures if a motif was not found.
                else:
                    # Subracts one from large_space_count if present, since
                    # a motif was not found.
                    if large_space_count > 0:
                        large_space_count -= 1  
                    # Subtracts one from double motif count if
                    # a motif has been found previously in the last
                    # four nucleotides.
                    if double_motif_count > 0:
                        double_motif_count -= 1
                    # Changes the last nucleotides if the spacing is higher
                    # than six nucleotides, indicating a new motif group. 
                    if motif_spacing > 6 and len(motifs_group_string) > 0:
                        # Grabs last seven nucleotides to add them to list.
                        # This signifies the end of a group.
                        last_nucleotides = around_motif[-7:]
                        last_index = around_index[-7:]
                        for nucleotide in last_nucleotides:
                            motifs_group_string += nucleotide
                            nucleotide_string += nucleotide
                        for index in last_index:
                            structure_index.append(index)    
                        # Adds string to list.
                        motifs_group_list.append(motifs_group_string)
                        # Resets string for new group.
                        motifs_group_string = ""
                    # Count to grab past groups if motif found
                    # within seven nucleotides. This will remove
                    # the previous nucleotides just added if a new
                    # group was not actually present.
                    if motif_spacing == 7: 
                        large_space_count = 3
                    if motif_spacing < 7:
                        motif_counter += 1
            # If a motif was recently found, appends the last nucleotides
            # to the strings. This will occur when a motif is towards
            # the end of the sequence.
            if (0 < motif_spacing <= 7):
                last_nucleotides = around_motif[-motif_spacing:]
                last_index = around_index[-motif_spacing:]
                for nucleotide in last_nucleotides:
                    nucleotide_string += nucleotide
                    motifs_group_string += nucleotide
                for index in last_index:
                    structure_index.append(index)
            # Appends nucleotide strings to each list for each region.
            # Does not append empty strings since these would later be
            # counted as a group.    
            nucleotide_list.append(nucleotide_string) \
            if len(nucleotide_string) > 0 else None
            nucleotide_intra_list.append(nucleotide_intra_string) \
            if len(nucleotide_intra_string) > 0 else None
            nucleotide_inter_list.append(nucleotide_inter_string) \
            if len(nucleotide_inter_string) > 0 else None
            motifs_group_list.append(motifs_group_string) \
            if len(motifs_group_string) > 0 else None  
        elif line[0] == ">":
            region_name = line_clean.replace(">", "")
            # Resets line count
            line_count = 0
        else:
            # Loops through each index saved in index lists and 
            # retrieves the secondary structure prediction.
            for index in structure_index:
                structure.append(line_clean[index])
            for index in y_index:
                y_structure.append(line_clean[index])
            for index in gc_index:
                gc_structure.append(line_clean[index])
        region_dictionary[region_name] = [nucleotide_list, nucleotide_intra_list, 
                                          nucleotide_inter_list, motifs_group_list,
                                          region_size]
        structure_dictionary[region_name] = [structure, y_structure, gc_structure]
    return region_dictionary, motif_title, motif_list, structure_dictionary                   

def output_region_and_file_data(file_dictionary, output):
    """Outputs secondary structure data for each region and each file.

        Arguments:
        file_dictionary -- Dictionary of secondary structure predictions for
            motifs and surrounding nucleotides per input file. The data is a 
            dictionary with the region as the key and the list of secondary 
            structure predictions as the value, arranged in three lists.
            {file_name: {"region_name": [structure, y_structure, gc_structure]}}. 
        output -- Output directory.

        Output:    
        secondary_per_region_output -- CSV file of secondary structure data for
            each region.
            "Region_Name,Surrounding_Paired_Percent,Y_Paired_Percent,GC_Paired_Percent".
        secondary_per_file_output -- CSV file of secondary structure data for
            each sample.
            "File_Name,Surrounding_Paired_Percent,Y_Paired_Percent,GC_Paired_Percent".    
    """
    # Checks to see if output directory ends with a "/".
    if output[-1] != "/":
        output = output + "/"
    # Makes output file names.
    secondary_per_region_output = (output + "secondary_structure_per_region.csv")
    secondary_per_file_output = (output + "secondary_structure_per_file.csv")    
    # Opens output files.
    secondary_per_region_output_file = open(secondary_per_region_output, 'w')
    secondary_per_file_output_file = open(secondary_per_file_output, 'w')
    # Writes headers to output files.
    secondary_per_region_output_file.write("Region_Name,Surrounding_Paired_Percent,Y_Paired_Percent,GC_Paired_Percent\n")
    secondary_per_file_output_file.write("File_Name,Surrounding_Paired_Percent,Y_Paired_Percent,GC_Paired_Percent\n")
    # Set used to store regions so that only one is outputted per file.
    region_set = set()
    # Loops through each file in file_dictionary and outputs data.
    for file_name in file_dictionary:
        # Used to store secondary structure sums per file
        structure_paired_sum = 0
        y_paired_sum = 0
        gc_paired_sum = 0
        structured_total = 0
        y_total = 0
        gc_total = 0      
        # Loops through each region in file_dictionary and outputs data.
        for region_name in file_dictionary[file_name]:
            # Used to store sums per region.
            structure_paired_sum_region = 0
            y_paired_sum_region = 0
            gc_paired_sum_region = 0
            structured_total_region = 0
            y_total_region = 0
            gc_total_region = 0
            # Grabs lists of secondary structure data.
            structure_list = file_dictionary[file_name][region_name][0]
            y_structure_list = file_dictionary[file_name][region_name][1]
            gc_structure_list = file_dictionary[file_name][region_name][2]
            # Gets total of "(" and ")" of items in each list.
            structure_paired_sum_region += structure_list.count("(")
            structure_paired_sum_region += structure_list.count(")")
            structured_total_region = len(structure_list)
            gc_paired_sum_region += gc_structure_list.count("(")
            gc_paired_sum_region += gc_structure_list.count(")")
            gc_total_region = len(gc_structure_list)
            y_paired_sum_region += y_structure_list.count("(")
            y_paired_sum_region += y_structure_list.count(")")
            y_total_region = len(y_structure_list)
            # Outputs the percentage of paired nucleotides for each region.
            # Only outputs if region has not been output yet.
            if region_name not in region_set:
                secondary_per_region_output_file.write(region_name + ",")
                if structured_total_region != 0:
                    secondary_per_region_output_file.write(str(100 *(structure_paired_sum_region / structured_total_region)) + ",")
                else:
                    secondary_per_region_output_file.write("0,")
                if y_total_region != 0:
                    secondary_per_region_output_file.write(str(100 *(y_paired_sum_region / y_total_region)) + ",")
                else:
                    secondary_per_region_output_file.write("0,")          
                if gc_total_region != 0:
                    secondary_per_region_output_file.write(str(100 *(gc_paired_sum_region / gc_total_region)) + "\n")
                else:
                    secondary_per_region_output_file.write("0\n")    
                region_set.add(region_name)
            # Adds to sums for each file.
            structure_paired_sum += structure_paired_sum_region
            y_paired_sum += y_paired_sum_region
            gc_paired_sum += gc_paired_sum_region
            structured_total += structured_total_region
            y_total += y_total_region
            gc_total += gc_total_region            
        # Outputs the percentage of paired nucleotides for each file.
        secondary_per_file_output_file.write(file_name + ",")
        if structured_total != 0:
            secondary_per_file_output_file.write(str(100 *(structure_paired_sum / structured_total)) + ",")
        else:
            secondary_per_file_output_file.write("0,")
        if y_total != 0:
            secondary_per_file_output_file.write(str(100 *(y_paired_sum / y_total)) + ",")
        else:
            secondary_per_file_output_file.write("0,")
        if gc_total != 0:
            secondary_per_file_output_file.write(str(100 *(gc_paired_sum / gc_total)) + "\n")        
        else:
            secondary_per_file_output_file.write("0\n")

def init_argparse():
    """Initiates the use of argparse. Returns parsed
        command line arguments.

        Arguments:
        See command line argument descriptions below.

        Output:
        Returns parser object. Object is used for function input.
    """
    # Script description.
    parser = argparse.ArgumentParser(description = "Gets structure secondary motif statistics for input FASTA files.")        
    # Command line arguments and descriptions.
    parser.add_argument("-fl", "--fasta", action = "store", type = str, nargs='+', 
                        help="List of FASTA files of RBP binding regions with secondary \
                        structure predictions from ViennaRNA.", 
                        required = True)
    parser.add_argument("-k", "--kmer", action = "store", type = int, 
                        help="Size of kmers to be analyzed.",
                        required = True)
    parser.add_argument("-o", "--output", action = "store", type = str,
                        help="Prefix for output files.",
                        required = True)
    parser.add_argument("-m", "--motif", action = "store", type = str, nargs='+',
                        default=["YGCY"],
                        help="List of motifs to be analyzed. Default is YGCY.",
                        required = False)
    parser.add_argument("-t", "--title", action = "store", type = str,
                        default="YGCY",
                        help="Title for motif output file. Default is YGCY.",
                        required = False)                                                                              
    return parser

def main():
    """Main script function. Executes other functions.

        Arguments:
        region_fasta_list -- (-fl) List of FASTA files of RBP 
            binding regions with secondary structure 
            predictions from ViennaRNA.
        kmer_size -- (-k) Size of kmer to be analyzed.
        output_file_prefix -- (-o) Prefix for output files.
        motif_list -- (-m) List of motifs to be analyzed. Default is YGCY.
            (Optional).
        motif_title -- (-t) Title for motif output file. Default is YGCY.
            (Optional).
            
        Output:    
        secondary_per_region_output -- CSV file of secondary structure data for
            each region.
            "Region_Name,Surrounding_Paired_Percent,Y_Paired_Percent,GC_Paired_Percent".
        secondary_per_file_output -- CSV file of secondary structure data for
            each sample.
            "File_Name,Surrounding_Paired_Percent,Y_Paired_Percent,GC_Paired_Percent".
    """
    # Uses argparse to organize command line arguments.
    args_to_parse = init_argparse()
    parsed = args_to_parse.parse_args()
    # Dictionary used to save file information.
    file_dictionary = {}
    for fasta_file in parsed.fasta:
        # Grabs fasta_file name.
        if "/" in fasta_file:
            fasta_file_name = fasta_file.split("/")[-1]
        else:
            fasta_file_name = fasta_file    
        # Organizes nucleotides and motifs.
        region_dictionary, motif_title, motif_list, structure_dictionary = \
        organize_nucleotides(fasta_file, parsed.kmer, parsed.motif, parsed.title)
        file_dictionary[fasta_file] = structure_dictionary
    output_region_and_file_data(file_dictionary, parsed.output)    

# Executes main function.
if __name__ == "__main__":
    main()