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

    """
    # Used to store nucleotide data. See description above.
    region_dictionary = {}
    # Checks motif_list for "YGCY", the default. If found, replaces
    # motif_list populates with all YGCY variants. If not found,
    # the input motif list is used. Uses motif_title to name output.
    if motif_list == ["YGCY"]:
         motif_list = ["CGCC", "CGCT", "TGCC", "TGCT"]
         motif_title = "YGCY"
    else:
         motif_list = motif_list
         motif_title = motif_title
    # Used to store name of binding region.
    region_name = ""
    open_file = open(fasta_file_name, 'r')
    # Loops through each line of the FASTA file. If title line
    # title is saved, if nucleotides, motif counts are generated.        
    for line in open_file:
        line_clean = line.strip("\n")
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
        # Used to detect title line.
        if line[0] != ">":
            # Grabs region size.
            region_size = len(line_clean)
            # Nucleotide list. List used to keep track of nucleotides around
            # motif of interest. Does not contain motifs.
            around_motif = []
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
                # Adds nucleotides to around_motif list. If list is empty,
                # adds all nucleotides. If not, adds only the last nucleotide.
                if (motif_count == 0 and len(around_motif) == 0):
                    for nucleotide in kmer:
                        around_motif.append(nucleotide)
                    motif_spacing += 4
                else:
                    around_motif.append(kmer[-1])
                    motif_spacing += 1
                # If kmer is in motif_list, saves up to the last seven nucleotides
                # and appends to nucleotide_list. Also determines what other
                # nucleotides need to be appended to other lists.
                if kmer in motif_list:
                    # Removes last three nucleotides from around_motif list if
                    # there was a double motif (i.e., YGCYGCY). If not, removes
                    # the standard four nucleotides. Set double_motif_count_flag
                    # to true so only three nucleotides are added to motif_group_string.
                    if double_motif_count != 0:
                        del around_motif[-3:]
                        motif_spacing -= 3
                        double_motif_count_flag = True 
                    else:
                        del around_motif[-4:]
                        motif_spacing -= 4    
                    # Used for first motif in binding region.
                    if motif_count == 0:
                        # Grabs previous seven non-motif nucleotides
                        # and adds them to nucleotide strings.
                        last_nucleotides = around_motif[-7:]
                        for nucleotide in last_nucleotides:
                            nucleotide_string += nucleotide
                            motifs_group_string += nucleotide
                    # If previous motif has been found, checks for the 
                    # double_motif_count_flag. If true no nucleotides
                    # are added since they would have already been added.        
                    elif double_motif_count_flag != True:
                        # Checks to see if motif is within seven nucleotides,
                        # this means that the motif is in the same group. If
                        # not, adds nucleotides to nucleotide_inter_list.
                        last_nucleotides = around_motif[-motif_spacing:]
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
                            # Adds to nucleotides to applicable strings.
                            for nucleotide in last_nucleotides:
                                nucleotide_string += nucleotide
                                motifs_group_string += nucleotide
                                nucleotide_intra_string += nucleotide
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
                                # Resets count because a motif was found.
                                large_space_count = 0   
                            # Now adds the last seven nucleotieds to the
                            # appropriate lists.    
                            last_nucleotides = around_motif[-7:]
                            for nucleotide in last_nucleotides:
                                nucleotide_string += nucleotide
                                motifs_group_string += nucleotide                                    
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
                        for nucleotide in last_nucleotides:
                            motifs_group_string += nucleotide
                            nucleotide_string += nucleotide
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
                for nucleotide in last_nucleotides:
                    nucleotide_string += nucleotide
                    motifs_group_string += nucleotide
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
        else:
            region_name = line_clean.replace(">", "")
        region_dictionary[region_name] = [nucleotide_list, nucleotide_intra_list, 
                                          nucleotide_inter_list, motifs_group_list,
                                          region_size]
    return region_dictionary, motif_title, motif_list                      

def get_motif_counts(region_dictionary, motif_list, fasta_file_name):
    """Analyzes nucleotide dictionaries for motif and secondary nucleotide
        characteristics.
        
        Arguments:
        region_dictionary -- Dictionary of nucleotide lists for the nucleotides
            surrounding each motif and the motifs themselves. The data is arranged
            is a list of four lists and an integer. The descriptions of each data
            structure is listed below. See organize_nucleotides for details.
            "nucleotide_list"- List, string list of nucleotides, seven nucleotides
                around all motifs.
            "nucleotide_intra_list"- List, string list of nucleotides between motifs
                in a motif group. A motif group is defined as motifs within seven 
                nucleotides of each other.
            "nucleotide_inter_list"- List, string list of nucleotides between motifs
                separated by more than seven nucleotides.
            "motifs_group_list"- List, string list of motifs and surrounding 
                nucleotides within a group.
            "region_size"- Integer, size of binding region (last variable).    
            {"region_name": [nucleotide_list, nucleotide_intra_list,
                             nucleotide_inter_list, motifs_group_list,
                             region_size]}.
        motif_list -- List of motifs to be analyzed. Default is YGCY.
        fasta_file_name -- String, name of FASTA file being analyzed.

        Output:
        region_stats_dictionary -- Dictionary with motif and secondary
            statistics for each region. Dictionary is organized by region.
            See nucleotide_content and motif_content for details.
            {"region_name": [nucleotide_list, nucleotide_intra_list,
                             nucleotide_inter_list, motifs_group_list,
                             region_size]}. 
        file_stats_dictionary -- Dictionary with motif and secondary
            statistics for each input file. See aggregate_file_stats
            for details.
            {"file_name": [nucleotide_count, nucleotide_intra_count,
                           nucleotide_inter_count, motif_region_counts]}.      
    """
    # Used to store statistics for each region. See below for list details.
    # region_stats_dictionary[region] = [region_stats_list]
    region_stats_dictionary = {}
    # Iterates through each region in region_dictionary and calculates
    # motif and secondary nucleotide statistics.
    for region in region_dictionary:
        # Grabs lists from region_dictionary.
        nucleotide_list = region_dictionary[region][0]
        nucleotide_intra_list = region_dictionary[region][1]
        nucleotide_inter_list = region_dictionary[region][2]
        motifs_group_list = region_dictionary[region][3]
        region_size = region_dictionary[region][4]
        # Grabs nucleotide counts for nucleotide lists. Count objects are all
        # dictionaries. See functions for details.
        nucleotide_count = nucleotide_content(nucleotide_list)
        nucleotide_intra_count = nucleotide_content(nucleotide_intra_list)
        nucleotide_inter_count = nucleotide_content(nucleotide_inter_list)
        # Grabs motif count and nucleotide count for motifs_group_list.
        motif_region_counts = motif_content(motifs_group_list, motif_list, region_size)
        region_stats_dictionary[region] = [nucleotide_count, nucleotide_intra_count,
                                           nucleotide_inter_count, motif_region_counts]
    # Gets aggregate statistics for each file.    
    file_stats_dictionary = aggregate_file_stats(region_stats_dictionary, fasta_file_name)
    return region_stats_dictionary, file_stats_dictionary   

def nucleotide_content(nucleotide_list):
    """Calculates nucleotide content for each nucleotide in a list.

        Arguments:
        nucleotide_list -- List of nucleotide strings to be analyzed.

        Output:
        nucleotide_count -- Dictionary of nucleotide counts for each nucleotide
            in nucleotide_list. Dictionary is organized by nucleotide. Contains
            nucleotide count total, nucleotide percentage, and average nucleotide
            length.
            nucleotide_count = {"A": 2, "C": 2, "G": 4, "T": 2,
                                "A_percentage": 0.0, "C_percentage": 0.0, 
                                "G_percentage": 0.0, "T_percentage": 0.0,  
                                "Total": 0, "Average_length": 0.0}.
    """
    # Used to store nucleotide counts.
    nucleotide_count = {}
    # Used to grab length of each nucleotide string.
    nucleotide_length = []
    # Initializes counts for each nucleotide.
    for nucleotide in ["A", "C", "G", "T", "Total"]:
        nucleotide_count[nucleotide] = 0
    # Iterates through each nucleotide set in nucleotide_list and adds 
    # counts to nucleotide_count.
    for nucleotide_string in nucleotide_list:
        nucleotide_length.append(len(nucleotide_string))
        for nucleotide in nucleotide_string:
            nucleotide_count["Total"] += 1
            nucleotide_count[nucleotide] += 1
    # Calculates percentage for each nucleotide.
    for nucleotide in ["A", "C", "G", "T"]:
        nucleotide_count[nucleotide + "_percentage"] = (100 * (nucleotide_count[nucleotide] 
                                                        / nucleotide_count["Total"])) \
        if nucleotide_count["Total"] != 0 else 0     
    # Calculates average nucleotide length and adds to nucleotide_count.        
    nucleotide_average = (sum(nucleotide_length) / len(nucleotide_length)) \
    if len(nucleotide_length) != 0 else 0
    nucleotide_count["Average_length"] = nucleotide_average
    return nucleotide_count

def motif_content(motifs_group_list, motif_list, region_size):
    """First sums motif and nucleotide content and adds
        to a list, using aggregate_motif_counts the counts
        are aggregated and values calculated.
    
        Arguments:
        motifs_group_list -- List, string list of motifs and surrounding 
            nucleotides within a group.
        motif_list -- String List of motifs to be analyzed. Default is YGCY.
        region_size -- Integer, size of binding region.

        Output:
        motif_region_counts -- Dictionary of motif statistics for each region.
                See aggregate_motif_counts for details.

    """
    # Output dictionary. See Output description above.
    motif_region_counts = {}
    # List used to store motif_counts and motif_counts
    # in dictionaries for each motif group.
    # [{"CGCC": 0, "CGCT": 0, "TGCC": 0, 
    #   "TGCT": 0, "Total_motif": 0,
    #    "A": 0, "C": 0, "G": 0, "T": 0, 
    #    "Total_nucleotide": 0}]
    motif_count_list = []
    # Grabs kmer size from first item in motif_list.
    kmer_size = len(motif_list[0])
    # Count used to count the number of groups.
    group_count = 0
    # Iterates through each motif in motifs_group_list and adds 
    # counts to motif_counts and motif_counts.
    for group in motifs_group_list:
        # Used to store motif and nucleotide counts of the
        # surrounding nucleotides for whole binding region.
        # See below for contents.
        motif_counts = {}
        # List used to store last four nucleotides added.
        # Used to remove last four nucleotides if it is
        # part of a motif.
        last_four_nucleotides = []
        # Initializes counts for each motif and nucleotide.
        for item in ["CGCC", "CGCT", "TGCC", "TGCT", "Total_motif",
                      "A", "C", "G", "T", "Total_nucleotide"]:
            motif_counts[item] = 0
        group_count += 1
        for start in range(len(group) - kmer_size + 1):     
            kmer = group[start:start + kmer_size]
            # Used for first kmer in group. Iterates through
            # and adds the entire kmer to last_three_nucleotides
            # and motif_counts. If has been added,
            # only adds the last nucleotide.
            if len(last_four_nucleotides) == 0:
                for nucleotide in kmer:
                    last_four_nucleotides.append(nucleotide)
                    motif_counts[nucleotide] += 1
                    motif_counts["Total_nucleotide"] += 1
            else:
                last_four_nucleotides.append(kmer[-1])
                motif_counts[kmer[-1]] += 1
                motif_counts["Total_nucleotide"] += 1
            # If kmer is in motif_list, adds counts to motif_counts
            # and motif_counts. Also removes last four nucleotides
            # from the counts since they are part of the motif.  
            if kmer in motif_list:
                motif_counts["Total_motif"] += 1
                motif_counts[kmer] += 1
                # Removes last four nucleotides counts. 
                four_nucleotides = last_four_nucleotides[-4:]
                for nucleotide in four_nucleotides:
                    motif_counts[nucleotide] -= 1
                    motif_counts["Total_nucleotide"] -= 1    
        # Append motif_counts and motif_counts to lists.
        motif_count_list.append(motif_counts)
    # Aggregates motif_counts and motif_counts per group
    # for each region.
    motif_region_counts = aggregate_motif_counts(motif_count_list, region_size)    
    return motif_region_counts

def aggregate_motif_counts(motif_count_list, region_size):
    """Aggregates motif and nucleotide counts for each motif group
        and calculates needed statistics.
    
            Arguments:
            motif_count_list -- List of dictionaries of motif counts for each
                motif group in a region. Dictionary is organized by motif group.
                [{"CGCC": 0, "CGCT": 0, "TGCC": 0, 
                  "TGCT": 0, "Total_motif": 0,
                  "A": 0, "C": 0, "G": 0, "T": 0, 
                  "Total_nucleotide": 0}].
            region_size -- Integer, size of binding region.  
    
            Output:
            motif_region_counts -- Dictionary of motif statistics for each region.
                See below for details of entries.
                nucleotide + "_rate_per_motif_per_group" -- Nucleotide rate 
                    (of total nucleotides) of each nucleotide per motif per group number
                nucleotide + "_per_group" -- Nucleotide type per group.    
                nucleotide + "_per_" + motif -- Nucleotide type per motif type.
                motif + "_per_group" -- Motif fraction type per group.
                motif + "_per_region" -- Motif fraction type per region.
                "motif_per_group" -- Average motifs per group.
                "nucleotide_per_group" -- Average nucleotides per group.
                "Total_motifs_per_region_length" -- Total motifs per region length.
                "Total_nucleotides_per_region" -- Total surrounding nucleotides per region length.
                {A_per_motif_per_group: 0, C_per_motif_per_group: 0, G_per_motif_per_group: 0,
                 T_per_motif_per_group: 0, A_per_group: 0, C_per_group: 0, G_per_group: 0,
                 T_per_group: 0,
                 "A_per_CGCC": 0, "A_per_CGCT": 0, "A_per_TGCC": 0, "A_per_TGCT": 0,
                 "C_per_CGCC": 0, "C_per_CGCT": 0, "C_per_TGCC": 0, "C_per_TGCT": 0,
                 "G_per_CGCC": 0, "G_per_CGCT": 0, "G_per_TGCC": 0, "G_per_TGCT": 0,
                 "T_per_CGCC": 0, "T_per_CGCT": 0, "T_per_TGCC": 0, "T_per_TGCT": 0,
                 "CGCC_per_group": 0, "CGCT_per_group": 0, "TGCC_per_group": 0, "TGCT_per_group": 0,
                 "CGCC_per_region": 0, "CGCT_per_region": 0, "TGCC_per_region": 0, 
                 "TGCT_per_region": 0, "motif_per_group": 0, "nucleotide_per_group": 0,
                 "Total_motifs_per_region_length": 0, "Total_nucleotides_per_region_length": 0}.        
    """
    # Initializes counts for each motif and nucleotide.
    # motif_counts_loop is used for iterating through
    # items in the dicionary.
    # Dictionary is formatted as follows:
    # {"CGCC": 0, "CGCT": 0, "TGCC": 0, 
    #  "TGCT": 0, "Total_motif": 0,
    #  "A": 0, "C": 0, "G": 0, "T": 0, 
    #  "Total_nucleotide": 0}
    motif_counts = {}
    motif_counts_loop = {}
    for item in ["CGCC", "CGCT", "TGCC", "TGCT", "Total_motif",
                 "A", "C", "G", "T", "Total_nucleotide"]:
        motif_counts[item] = 0
        motif_counts_loop[item] = 0
    # Iterates through each motif group in motif_count_list and
    # motif_nucleotide_count_list and aggregates counts.
    for group_dictionary in motif_count_list:
        # Used to store motif and nucleotide counts of the
        # surrounding nucleotides for whole binding region.
        # Iterates through each item in group and adds 
        # counts to motif_counts.
        for item in motif_counts_loop:
            motif_counts[item] += group_dictionary[item]
    # Calculates nucleotide rate (of total nucleotides) of each nucleotide per motif per group number and
    # nucleotides per group.
    for nucleotide in ["A", "C", "G", "T"]:
        motif_counts[nucleotide + "_rate_per_motif_per_group"] = ((motif_counts[nucleotide] 
                                                                   / motif_counts["Total_nucleotide"])
                                                                    / (motif_counts["Total_motif"]
                                                                     / len(motif_count_list))) \
        if ((motif_counts["Total_nucleotide"] != 0) and (motif_counts["Total_motif"] != 0) 
            and (len(motif_count_list)) != 0) else 0
        # Separate statement.                                                                
        motif_counts[nucleotide + "_per_group"] = (motif_counts[nucleotide] 
                                                   / len(motif_count_list)) \
        if len(motif_count_list) != 0 else 0   
    # Calculates each nucleotide type per motif type.
    for motif in ["CGCC", "CGCT", "TGCC", "TGCT"]:
        for nucleotide in ["A", "C", "G", "T"]:
            motif_counts[nucleotide + "_per_" + motif] = (motif_counts[nucleotide] 
                                                          / motif_counts[motif]) \
            if motif_counts[motif] != 0 else 0 
    # Calculates motif fraction type per group and motif fraction type per region.
    for motif in ["CGCC", "CGCT", "TGCC", "TGCT"]:
        motif_counts[motif + "_per_group"] = (motif_counts[motif] / len(motif_count_list)) \
        if len(motif_count_list) != 0 else 0 
        # Separate statement.  
        motif_counts[motif + "_per_region"] = (motif_counts[motif] 
                                               / motif_counts["Total_motif"]) \
        if motif_counts["Total_motif"] != 0 else 0                                          
    # Calculates total motifs and nucleotides per group.
    if len(motif_count_list) != 0:
        motif_counts["motif_per_group"] = (motif_counts["Total_motif"] / len(motif_count_list))
        motif_counts["nucleotide_per_group"] = (motif_counts["Total_nucleotide"] / len(motif_count_list))
    else:
        motif_counts["motif_per_group"] = 0
        motif_counts["nucleotide_per_group"] = 0
    # Saves total motifs per region and surrounding nucleotides per region size.
    if region_size != 0:
        motif_counts["Total_motifs_per_region_length"] = 100 *(motif_counts["Total_motif"] / region_size)
        motif_counts["Total_nucleotides_per_region_length"] = 100 *(motif_counts["Total_nucleotide"] / region_size)
    else:
        motif_counts["Total_motifs_per_region_length"] = 0
        motif_counts["Total_nucleotides_per_region_length"] = 0         
    return motif_counts

def aggregate_file_stats(region_stats_dictionary, fasta_file_name):
    """Aggregates motif and secondary nucleotide statistics for each file.
    
        Arguments:
        region_stats_dictionary -- Dictionary with motif and seconardary
            statistics for each region. Dictionary is organized by region.
            For nucleotide_count, nucleotide_intra_count, nucleotide_inter_count,
            see nucleotide_content for details, for motif_region_counts, see
            aggregate_motif_counts for details.
            {region: [nucleotide_count, nucleotide_intra_count,
                      nucleotide_inter_count, motif_region_counts]}. 
        fasta_file_name -- String, name of FASTA file being analyzed.                             
    
        Output:
        file_stats_dictionary -- Dictionary with motif and secondary
            statistics for each input file. Same format as region_stats_dictionary.
            but contains statistics for entire files instead of regions.     
    """
    # Used to store statistics for each file. See below for list details.
    # file_stats_dictionary[file_name] = [file_stats_list]
    file_stats_dictionary = {}
    # Dictionaries used to store aggregate statistics for each file.
    nucleotide_count_file = {}
    nucleotide_intra_count_file = {}
    nucleotide_inter_count_file = {}
    motif_region_counts_file = {}
    # Count used to count the number of regions.
    region_count = 0
    # Iterates through each region in region_stats_dictionary and calculates
    # motif and secondary nucleotide statistics.
    for region in region_stats_dictionary:
        region_count += 1 
        # Grabs lists from region_stats_dictionary.
        nucleotide_count = region_stats_dictionary[region][0]
        nucleotide_intra_count = region_stats_dictionary[region][1]
        nucleotide_inter_count = region_stats_dictionary[region][2]
        motif_region_counts = region_stats_dictionary[region][3]
        # Calculates percentage sum for each nucleotide list. Nucleotide lists
        # and motif list calculations separated for readability.
        for nucleotide in ["A", "C", "G", "T"]:
            if (nucleotide + "_percentage_sum") not in nucleotide_count_file:
                nucleotide_count_file[nucleotide + "_percentage_sum"] = nucleotide_count[nucleotide + "_percentage"]
                nucleotide_intra_count_file[nucleotide + "_percentage_sum"] = nucleotide_intra_count[nucleotide + "_percentage"]
                nucleotide_inter_count_file[nucleotide + "_percentage_sum"] = nucleotide_inter_count[nucleotide + "_percentage"]
            else:
                nucleotide_count_file[nucleotide + "_percentage_sum"] += nucleotide_count[nucleotide + "_percentage"]
                nucleotide_intra_count_file[nucleotide + "_percentage_sum"] += nucleotide_intra_count[nucleotide + "_percentage"]
                nucleotide_inter_count_file[nucleotide + "_percentage_sum"] += nucleotide_inter_count[nucleotide + "_percentage"]
        # Uses average nucleotide length and adds to nucleotide_count to get sum per file.
        if ("Average_length_sum") not in nucleotide_count_file:
            nucleotide_count_file["Average_length_sum"] = nucleotide_count["Average_length"]
            nucleotide_intra_count_file["Average_length_sum"] = nucleotide_intra_count["Average_length"]
            nucleotide_inter_count_file["Average_length_sum"] = nucleotide_inter_count["Average_length"]
        else:
            nucleotide_count_file["Average_length_sum"] += nucleotide_count["Average_length"]
            nucleotide_intra_count_file["Average_length_sum"] += nucleotide_intra_count["Average_length"]
            nucleotide_inter_count_file["Average_length_sum"] += nucleotide_inter_count["Average_length"]
        # Iterates through required motif statistics in motif_region_counts, sums all the counts, and then
        # calculates the average.
        # Calculates sum of nucleotide rate (of total nucleotides) per motif per group number, per group, and
        # nucleotide rate per region.
        for nucleotide in ["A", "C", "G", "T"]:
            if ((nucleotide + "_per_motif_per_group" + "_sum") not in motif_region_counts_file and 
                (nucleotide + "_per_group_sum") not in motif_region_counts_file):
                motif_region_counts_file[nucleotide + "_rate_per_motif_per_group_sum"] = motif_region_counts[nucleotide + "_rate_per_motif_per_group"]
                motif_region_counts_file[nucleotide + "_per_group_sum"] = motif_region_counts[nucleotide + "_per_group"]
            else:
                motif_region_counts_file[nucleotide + "_rate_per_motif_per_group_sum"] += motif_region_counts[nucleotide + "_rate_per_motif_per_group"]
                motif_region_counts_file[nucleotide + "_per_group_sum"] += motif_region_counts[nucleotide + "_per_group"]        
        # Calculates the sum of each nucleotide type per motif type.
        for motif in ["CGCC", "CGCT", "TGCC", "TGCT"]:
            for nucleotide in ["A", "C", "G", "T"]:
                if (nucleotide + "_per_" + motif + "_sum") not in motif_region_counts_file:
                    motif_region_counts_file[nucleotide + "_per_" + motif + "_sum"] = motif_region_counts[nucleotide + "_per_" + motif]
                else:
                    motif_region_counts_file[nucleotide + "_per_" + motif + "_sum"] += motif_region_counts[nucleotide + "_per_" + motif]
        # Calculates sum of motifs per group and motifs per region.
        for motif in ["CGCC", "CGCT", "TGCC", "TGCT"]:
            if ((motif + "_per_group" + "_sum") not in motif_region_counts_file and
                (motif + "_per_region" + "_sum") not in motif_region_counts_file):
                motif_region_counts_file[motif + "_per_group_sum"] = motif_region_counts[motif + "_per_group"]
                motif_region_counts_file[motif + "_per_region_sum"] = motif_region_counts[motif + "_per_region"]
            else:
                motif_region_counts_file[motif + "_per_group_sum"] += motif_region_counts[motif + "_per_group"]
                motif_region_counts_file[motif + "_per_region_sum"] += motif_region_counts[motif + "_per_region"]       
        # Calculates sum of total motifs, nucleotides per group, sum of total motifs, and nucleotides per region.
        if (("motif_per_group_sum") not in motif_region_counts_file and
            ("nucleotide_per_group_sum") not in motif_region_counts_file and
            ("Total_motifs_per_region_length_sum") not in motif_region_counts_file and
            ("Total_nucleotides_per_region_length_sum") not in motif_region_counts_file):
            motif_region_counts_file["motif_per_group_sum"] = motif_region_counts["motif_per_group"]
            motif_region_counts_file["nucleotide_per_group_sum"] = motif_region_counts["nucleotide_per_group"]
            motif_region_counts_file["Total_motifs_per_region_length_sum"] = motif_region_counts["Total_motifs_per_region_length"]
            motif_region_counts_file["Total_nucleotides_per_region_length_sum"] = motif_region_counts["Total_nucleotides_per_region_length"]
        else:
            motif_region_counts_file["motif_per_group_sum"] += motif_region_counts["motif_per_group"]
            motif_region_counts_file["nucleotide_per_group_sum"] += motif_region_counts["nucleotide_per_group"]
            motif_region_counts_file["Total_motifs_per_region_length_sum"] += motif_region_counts["Total_motifs_per_region_length"]
            motif_region_counts_file["Total_nucleotides_per_region_length_sum"] += motif_region_counts["Total_nucleotides_per_region_length"]      
    # Calculates averages of all the above summed statistics.
    # Nucleotide files.
    if region_count != 0:
        for nucleotide in ["A", "C", "G", "T"]:
            nucleotide_count_file[nucleotide + "_percentage"] = (nucleotide_count_file[nucleotide + "_percentage_sum"] 
                                                                 / region_count)
            nucleotide_intra_count_file[nucleotide + "_percentage"] = (nucleotide_intra_count_file[nucleotide + "_percentage_sum"] 
                                                                       / region_count)
            nucleotide_inter_count_file[nucleotide + "_percentage"] = (nucleotide_inter_count_file[nucleotide + "_percentage_sum"] 
                                                                       / region_count)
        nucleotide_count_file["Average_length"] = (nucleotide_count_file["Average_length_sum"]
                                                   / region_count)
        nucleotide_intra_count_file["Average_length"] = (nucleotide_intra_count_file["Average_length_sum"]
                                                         / region_count)
        nucleotide_inter_count_file["Average_length"] = (nucleotide_inter_count_file["Average_length_sum"]
                                                         / region_count)
        # Motif files.
        for nucleotide in ["A", "C", "G", "T"]:
            motif_region_counts_file[nucleotide + "_rate_per_motif_per_group"] = (motif_region_counts_file[nucleotide + "_rate_per_motif_per_group_sum"]
                                                                                  / region_count)
            motif_region_counts_file[nucleotide + "_per_group"] = (motif_region_counts_file[nucleotide + "_per_group_sum"]
                                                                   / region_count)
        for motif in ["CGCC", "CGCT", "TGCC", "TGCT"]:
            for nucleotide in ["A", "C", "G", "T"]:
                motif_region_counts_file[nucleotide + "_per_" + motif] = (motif_region_counts_file[nucleotide + "_per_" + motif + "_sum"] 
                                                                          / region_count)
        for motif in ["CGCC", "CGCT", "TGCC", "TGCT"]:
            motif_region_counts_file[motif + "_per_group"] = (motif_region_counts_file[motif + "_per_group_sum"] 
                                                              / region_count)
            motif_region_counts_file[motif + "_per_region"] = (motif_region_counts_file[motif + "_per_region_sum"] 
                                                               / region_count)
        motif_region_counts_file["motif_per_group"] = (motif_region_counts_file["motif_per_group_sum"]  
                                                       / region_count)
        motif_region_counts_file["nucleotide_per_group"] = (motif_region_counts_file["nucleotide_per_group_sum"]
                                                            / region_count)
        motif_region_counts_file["Total_motifs_per_region_length"] = (motif_region_counts_file["Total_motifs_per_region_length_sum"]
                                                                      / region_count)
        motif_region_counts_file["Total_nucleotides_per_region_length"] = (motif_region_counts_file["Total_nucleotides_per_region_length_sum"]
                                                                           / region_count)  
    # Adds each dictionary to a list and saves in file_stats_dictionary.
    file_stats_list = [nucleotide_count_file, nucleotide_intra_count_file,
                       nucleotide_inter_count_file, motif_region_counts_file]
    # Loops through each dictionary and removes all keys with the phrase "sum".
    for dictionary in file_stats_list:
        for key in list(dictionary.keys()):
            if "sum" in key:
                del dictionary[key] 
    file_stats_dictionary[fasta_file_name] = file_stats_list
    return file_stats_dictionary


def organize_files(region_fasta_list, kmer_size, motif_list, motif_title,
                   output_file_prefix):
    """Organizes input files and executes them through the required
        functions.

        Arguments:
        region_fasta_list -- List of FASTA files of RBP binding regions.
            Should include previously generated random regions from
            make_random_bed_for_motif.py. The random files are detected
            by searching for "random" in the file name.
        kmer_size -- Size of kmer to be analyzed.
        motif_list -- List of motifs to be analyzed. Default is YGCY.
        motif_title -- Title for motif output file. Default is YGCY.
        output_file_prefix -- Prefix for output files.

        Output:
        Uses output_stats_per_region and output_stats_per_file to 
        output statistics for each region and file.

    """
    # Used to save all file statistics when looping
    file_stats_dictionary_all = {}
    # Iterates through each FASTA file and executes functions.
    for fasta_file_name in region_fasta_list:
        # Grabs file name from region_fasta_file.
        if fasta_file_name.find("/") == -1:
            file_name = fasta_file_name
        else:    
            file_name = fasta_file_name.split("/")[-1]
        # Grabs intitial informaion for nucleotides and motifs.    
        (region_dictionary, motif_title, motif_list) = organize_nucleotides(fasta_file_name, kmer_size, 
                                                                            motif_list, motif_title)
        region_stats_dictionary, file_stats_dictionary = get_motif_counts(region_dictionary, motif_list, 
                                                                          fasta_file_name)
        # Adds file_stats_dictionary to file_stats_dictionary_all.
        file_stats_dictionary_all.update(file_stats_dictionary)
        # Outputs statistics for each region on a per file basis.
        output_stats_per_region(region_stats_dictionary, output_file_prefix, file_name)
    # Outputs all file statistics into a single output file.    
    output_stats_per_file(file_stats_dictionary_all, output_file_prefix) 

def output_stats_per_region(region_stats_dictionary, output_file_prefix, file_name):
    """Outputs statistics for each region.
    
            Arguments:
            region_stats_dictionary -- Dictionary with motif and seconardary
                statistics for each region. Dictionary is organized by region.
                For nucleotide_count, nucleotide_intra_count, nucleotide_inter_count,
                see nucleotide_content for details, for motif_region_counts, see
                aggregate_motif_counts for details.
                {region: [nucleotide_count, nucleotide_intra_count,
                        nucleotide_inter_count, motif_region_counts]}.
            output_file_prefix -- Prefix for output files.
            file_name -- Name of file being analyzed.
    
            Output:
            CSV files -- CSV file of nucleotide and motif secondary statistics.
            Contains all the statisitics in region_stats_dictionary. The titles
            of the file and an explanatation is listed below. The files are split
            into two files to make it easier to read.
            File 1- "_secondary_motif_stats_per_region_1.csv"
            Average percentages of each nucleotides surrounding motifs, within motif groups,
            and between motif groups and average length of surrounding nucleotides:
                All surrounding- All_A_%, All_C_%, All_G_%, All_T_%, 
                    All_average_length
                Within motif group- Intra_A_%, Intra_C_%, Intra_G_%, Intra_T_%, 
                    Intra_average_length
                Between motif group- Inter_A_%, Inter_C_%, Inter_G_%, Inter_T_%, 
                    Inter_average_length
            Average nucleotide rate of each nucleotide per motif per group number:
                A_rate_per_motif_per_group, C_rate_per_motif_per_group, 
                    G_rate_per_motif_per_group, T_rate_per_motif_per_group        
            Nucleotide types per group average:
                A_per_group, C_per_group, G_per_group, T_per_group
            File 2- "_secondary_motif_stats_per_region_2.csv"    
            Surrounding nucleotide type per motif type:
                A_per_CGCC, A_per_CGCT, A_per_TGCC, A_per_TGCT, C_per_CGCC, 
                    C_per_CGCT, C_per_TGCC, C_per_TGCT, G_per_CGCC, G_per_CGCT, 
                    G_per_TGCC, G_per_TGCT, T_per_CGCC, T_per_CGCT, T_per_TGCC, 
                    T_per_TGCT
            Motif types per group average:
                CGCC_per_group, CGCT_per_group, TGCC_per_group, TGCT_per_group
             Average motif fraction per region-
                CGCC_per_region, CGCT_per_region, TGCC_per_region, TGCT_per_region    
            Average motifs per group and average surrounding nucleotides per group
            (this is the number per 100 nucleotides):
                Motif_per_group, Nucleotide_per_group
            Total average motifs per region and total average surrounding 
                nucleotides per region length:
                Total_motifs_per_region_length, Total_nucleotides_per_region_length
        """
    # Creates output file names.
    output_file_name_1 = f"{output_file_prefix}_{file_name}_secondary_motif_stats_per_region_1.csv"
    output_file_name_2 = f"{output_file_prefix}_{file_name}_secondary_motif_stats_per_region_2.csv"
    # Creates output files.
    output_file_1 = open(output_file_name_1, "w")
    output_file_2 = open(output_file_name_2, "w")
    # Writes titles to output files.
    output_file_1.write("Region,All_A_%,All_C_%,All_G_%,All_T_%,All_average_length,"
                         "Intra_A_%,Intra_C_%,Intra_G_%,Intra_T_%,Intra_average_length,"
                         "Inter_A_%,Inter_C_%,Inter_G_%,Inter_T_%,Inter_average_length,"
                         "A_rate_per_motif_per_group,C_rate_per_motif_per_group,"
                         "G_rate_per_motif_per_group,T_rate_per_motif_per_group,"
                         "A_per_group,C_per_group,G_per_group,T_per_group\n")
    output_file_2.write("Region,A_per_CGCC,A_per_CGCT,A_per_TGCC,A_per_TGCT,C_per_CGCC,"
                         "C_per_CGCT,C_per_TGCC,C_per_TGCT,G_per_CGCC,G_per_CGCT,"
                         "G_per_TGCC,G_per_TGCT,T_per_CGCC,T_per_CGCT,T_per_TGCC,"
                         "T_per_TGCT,CGCC_per_group,CGCT_per_group,TGCC_per_group,TGCT_per_group,"
                         "CGCC_per_region,CGCT_per_region,TGCC_per_region,TGCT_per_region,"
                         "Motif_per_group,Nucleotide_per_group,"
                         "Total_motifs_per_region_length,Total_nucleotides_per_region_length\n")
    # Iterates through each region in region_stats_dictionary and outputs
    # statistics.
    for region in region_stats_dictionary:
        # Grabs lists from region_stats_dictionary.
        nucleotide_count = region_stats_dictionary[region][0]
        nucleotide_intra_count = region_stats_dictionary[region][1]
        nucleotide_inter_count = region_stats_dictionary[region][2]
        motif_region_counts = region_stats_dictionary[region][3]
        # Writes region name to output files.
        output_file_1.write(region + ",")
        output_file_2.write(region + ",")
        # Writes average percentages of each nucleotides surrounding motifs, within motif groups,
        # and between motif groups and average length of surrounding nucleotides.
        # All surrounding.
        for nucleotide in ["A", "C", "G", "T"]:
            output_file_1.write(str(nucleotide_count[nucleotide + "_percentage"]) + ",")
        output_file_1.write(str(nucleotide_count["Average_length"]) + ",")
        # Within motif group.
        for nucleotide in ["A", "C", "G", "T"]:
            output_file_1.write(str(nucleotide_intra_count[nucleotide + "_percentage"]) + ",")
        output_file_1.write(str(nucleotide_intra_count["Average_length"]) + ",")
        # Between motif group.
        for nucleotide in ["A", "C", "G", "T"]:
            output_file_1.write(str(nucleotide_inter_count[nucleotide + "_percentage"]) + ",")
        output_file_1.write(str(nucleotide_inter_count["Average_length"]) + ",")
        # Writes average nucleotide rate of each nucleotide per motif per group number.
        for nucleotide in ["A", "C", "G", "T"]:
            output_file_1.write(str(motif_region_counts[nucleotide + "_rate_per_motif_per_group"]) + ",")
        # Writes nucleotide types per group average.
        for nucleotide in ["A", "C", "G", "T"]:
            output_file_1.write(str(motif_region_counts[nucleotide + "_per_group"]) + ",")    
        output_file_1.write("\n")    
        # Writes surrounding nucleotide type per motif type.
        for nucleotide in ["A", "C", "G", "T"]:
            for motif in ["CGCC", "CGCT", "TGCC", "TGCT"]:
                output_file_2.write(str(motif_region_counts[nucleotide + "_per_" + motif]) + ",")
        # Writes motif types per group average.
        for motif in ["CGCC", "CGCT", "TGCC", "TGCT"]:
            output_file_2.write(str(motif_region_counts[motif + "_per_group"]) + ",")
        # Writes average motif fraction per region.
        for motif in ["CGCC", "CGCT", "TGCC", "TGCT"]:
            output_file_2.write(str(motif_region_counts[motif + "_per_region"]) + ",")
        # Writes average motifs per group and average surrounding nucleotides per group.
        output_file_2.write(str(motif_region_counts["motif_per_group"]) + ",")
        output_file_2.write(str(motif_region_counts["nucleotide_per_group"]) + ",")
        # Writes total average motifs per region and total average surrounding nucleotides per region length.
        output_file_2.write(str(motif_region_counts["Total_motifs_per_region_length"]) + ",")
        output_file_2.write(str(motif_region_counts["Total_nucleotides_per_region_length"]) + "\n")
    # Closes output files.
    output_file_1.close()
    output_file_2.close()

def output_stats_per_file(file_stats_dictionary_all, output_file_prefix):
    """Outputs statistics for each file.
        
                Arguments:
                file_stats_dictionary_all -- Dictionary with motif and seconardary
                    statistics for each input file.
                output_file_prefix -- Prefix for output files.
        
                Output:
                CSV files -- CSV file of nucleotide and motif secondary statistics.
                Contains all the statisitics in file_stats_dictionary_all. The titles
                of the file and an explanatation is listed below. The files are split
                into two files to make it easier to read.
                File 1- "_secondary_motif_stats_per_file_1.csv"
                Average percentages of each nucleotides surrounding motifs, within motif groups,
                and between motif groups and average length of surrounding nucleotides:
                    All_A_%, All_C_%, All_G_%, All_T_%, 
                        All_average_length
                    Within motif group- Intra_A_%, Intra_C_%, Intra_G_%, Intra_T_%, 
                        Intra_average_length
                    Between motif group- Inter_A_%, Inter_C_%, Inter_G_%, Inter_T_%, 
                        Inter_average_length
                Average nucleotide rate of each nucleotide per motif per group number:
                    A_rate_per_motif_per_group, C_rate_per_motif_per_group, 
                        G_rate_per_motif_per_group, T_rate_per_motif_per_group        
                Nucleotide types per group average:
                    A_per_group, C_per_group, G_per_group, T_per_group
                File 2- "_secondary_motif_stats_per_file_2.csv"    
                Surrounding nucleotide type per motif type:
                    A_per_CGCC, A_per_CGCT, A_per_TGCC, A_per_TGCT, C_per_CGCC, 
                        C_per_CGCT, C_per_TGCC, C_per_TGCT, G_per_CGCC, G_per_CGCT, 
                        G_per_TGCC, G_per_TGCT, T_per_CGCC, T_per_CGCT, T_per_TGCC, 
                        T_per_TGCT
                Motif types per group average:
                    CGCC_per_group, CGCT_per_group, TGCC_per_group, TGCT_per_group
                Average motif fraction per region:
                    CGCC_per_region, CGCT_per_region, TGCC_per_region, TGCT_per_region    
                Average motifs per group and average surrounding nucleotides per group
                (this is the number per 100 nucleotides):
                    Motif_per_group, Nucleotide_per_group
                Total average motifs per region and total average surrounding
                    nucleotides per region length:
                    Total_motifs_per_region_length, Total_nucleotides_per_region_length
                Total example dictionary below (values arbitrary):
                {'CGCC': 1, 'CGCT': 1, 'TGCC': 1, 'TGCT': 1, 'Total_motif': 4,
                 'A': 2, 'C': 1, 'G': 1, 'T': 6, 'Total_nucleotide': 10,
                 'A_rate_per_motif_per_group': 0.05, 'A_per_group': 2.0, 
                 'C_rate_per_motif_per_group': 0.025, 'C_per_group': 1.0, 
                 'G_rate_per_motif_per_group': 0.025, 'G_per_group': 1.0, 
                 'T_rate_per_motif_per_group': 0.15, 'T_per_group': 6.0, 
                 'A_per_CGCC': 2.0, 'C_per_CGCC': 1.0, 'G_per_CGCC': 1.0, 
                 'T_per_CGCC': 6.0, 'A_per_CGCT': 2.0, 'C_per_CGCT': 1.0, 
                 'G_per_CGCT': 1.0, 'T_per_CGCT': 6.0, 'A_per_TGCC': 2.0, 
                 'C_per_TGCC': 1.0, 'G_per_TGCC': 1.0, 'T_per_TGCC': 6.0, 
                 'A_per_TGCT': 2.0, 'C_per_TGCT': 1.0, 'G_per_TGCT': 1.0, 
                 'T_per_TGCT': 6.0, 'CGCC_per_group': 1.0, 'CGCC_per_region': 0.25, 
                 'CGCT_per_group': 1.0, 'CGCT_per_region': 0.25, 'TGCC_per_group': 1.0, 
                 'TGCC_per_region': 0.25, 'TGCT_per_group': 1.0, 'TGCT_per_region': 0.25, 
                 'motif_per_group': 4.0, 'nucleotide_per_group': 10.0, 
                 'Total_motifs_per_region_length': 4.0, 'Total_nucleotides_per_region_length': 10.0}    
    """
    # Creates output file names.
    output_file_name_1 = (output_file_prefix + "_secondary_motif_stats_per_file_1.csv")
    output_file_name_2 = (output_file_prefix + "_secondary_motif_stats_per_file_2.csv")
    # Creates output files.
    output_file_1 = open(output_file_name_1, "w")
    output_file_2 = open(output_file_name_2, "w")
    # Writes titles to output files.
    output_file_1.write("File_Name,All_A_%,All_C_%,All_G_%,All_T_%,All_average_length,"
                         "Intra_A_%,Intra_C_%,Intra_G_%,Intra_T_%,Intra_average_length,"
                         "Inter_A_%,Inter_C_%,Inter_G_%,Inter_T_%,Inter_average_length,"
                         "A_rate_per_motif_per_group,C_rate_per_motif_per_group,"
                         "G_rate_per_motif_per_group,T_rate_per_motif_per_group,"
                         "A_per_group,C_per_group,G_per_group,T_per_group\n")
    output_file_2.write("File_Name,A_per_CGCC,A_per_CGCT,A_per_TGCC,A_per_TGCT,C_per_CGCC,"
                         "C_per_CGCT,C_per_TGCC,C_per_TGCT,G_per_CGCC,G_per_CGCT,"
                         "G_per_TGCC,G_per_TGCT,T_per_CGCC,T_per_CGCT,T_per_TGCC,"
                         "T_per_TGCT,CGCC_per_group,CGCT_per_group,TGCC_per_group,TGCT_per_group,"
                         "CGCC_per_region,CGCT_per_region,TGCC_per_region,TGCT_per_region,"
                         "Motif_per_group,Nucleotide_per_group,"
                         "Total_motifs_per_region_length,Total_nucleotides_per_region_length\n")
    # Iterates through each file in file_stats_dictionary_all and outputs
    # statistics.
    for file_name in file_stats_dictionary_all:
        # Grabs lists from file_stats_dictionary_all.
        nucleotide_count_file = file_stats_dictionary_all[file_name][0]
        nucleotide_intra_count_file = file_stats_dictionary_all[file_name][1]
        nucleotide_inter_count_file = file_stats_dictionary_all[file_name][2]
        motif_region_counts_file = file_stats_dictionary_all[file_name][3]
        # Writes file name to output files.
        output_file_1.write(file_name + ",")
        output_file_2.write(file_name + ",")
        # Writes average percentages of each nucleotides surrounding motifs, within motif groups,
        # and between motif groups and average length of surrounding nucleotides.
        # All surrounding.
        for nucleotide in ["A", "C", "G", "T"]:
            output_file_1.write(str(nucleotide_count_file[nucleotide + "_percentage"]) + ",")
        output_file_1.write(str(nucleotide_count_file["Average_length"]) + ",")
        # Within motif group.
        for nucleotide in ["A", "C", "G", "T"]:
            output_file_1.write(str(nucleotide_intra_count_file[nucleotide + "_percentage"]) + ",")
        output_file_1.write(str(nucleotide_intra_count_file["Average_length"]) + ",")
        # Between motif group.
        for nucleotide in ["A", "C", "G", "T"]:
            output_file_1.write(str(nucleotide_inter_count_file[nucleotide + "_percentage"]) + ",")
        output_file_1.write(str(nucleotide_inter_count_file["Average_length"]) + ",")
        # Writes average nucleotide rate of each nucleotide per motif per group number.
        for nucleotide in ["A", "C", "G", "T"]:
            output_file_1.write(str(motif_region_counts_file[nucleotide + "_rate_per_motif_per_group"]) + ",")
        # Writes nucleotide types per group average.
        for nucleotide in ["A", "C", "G", "T"]:
            output_file_1.write(str(motif_region_counts_file[nucleotide + "_per_group"]) + ",")
        output_file_1.write("\n")    
        # Writes surrounding nucleotide type per motif type.
        for nucleotide in ["A", "C", "G", "T"]:
            for motif in ["CGCC", "CGCT", "TGCC", "TGCT"]:
                output_file_2.write(str(motif_region_counts_file[nucleotide + "_per_" + motif]) + ",")
        # Writes motif types per group average.
        for motif in ["CGCC", "CGCT", "TGCC", "TGCT"]:
            output_file_2.write(str(motif_region_counts_file[motif + "_per_group"]) + ",")
        # Writes average motif fraction per region.
        for motif in ["CGCC", "CGCT", "TGCC", "TGCT"]:
            output_file_2.write(str(motif_region_counts_file[motif + "_per_region"]) + ",")
        # Writes average motifs per group and average surrounding nucleotides per group.
        output_file_2.write(str(motif_region_counts_file["motif_per_group"]) + ",")
        output_file_2.write(str(motif_region_counts_file["nucleotide_per_group"]) + ",")
        # Writes total average motifs per region and total average surrounding nucleotides per region length.
        output_file_2.write(str(motif_region_counts_file["Total_motifs_per_region_length"]) + ",")
        output_file_2.write(str(motif_region_counts_file["Total_nucleotides_per_region_length"]) + "\n")
    # Closes output files.
    output_file_1.close()
    output_file_2.close()

def init_argparse():
    """Initiates the use of argparse. Returns parsed
        command line arguments.

        Arguments:
        See command line argument descriptions below.

        Output:
        Returns parser object. Object is used for function input.
    """
    # Script description.
    parser = argparse.ArgumentParser(description = "Gets secondary motif statistics for input FASTA files.")        
    # Command line arguments and descriptions.
    parser.add_argument("-fl", "--fasta", action = "store", type = str, nargs='+', 
                        help="List of FASTA files of RBP binding regions. Random files used \
                              for normalization should be included and contain random in the \
                              name.", 
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
        region_fasta_list -- (-fl) List of FASTA files of RBP binding 
            regions. Random files used for normalization should be 
            included and contain "random" in the name.
        kmer_size -- (-k) Size of kmer to be analyzed.
        output_file_prefix -- (-o) Prefix for output files.
        motif_list -- (-m) List of motifs to be analyzed. Default is YGCY.
            (Optional).
        motif_title -- (-t) Title for motif output file. Default is YGCY.
            (Optional).
            
        Output:    
        Output file CSV -- CSV file of nucleotide and motif secondary statistics.
            Contains all the statisitics in file_stats_dictionary_all. The titles
            of the file and an explanation is listed in output_stats_per_file. 
            The files are split into two files to make it easier to read. See 
            output_stats_per_file for details.
        Output region CSV -- CSV file of nucleotide and motif secondary statistics.
            Contains all the statisitics in region_stats_dictionary. The titles
            of the file and an explanation is listed in output_stats_per_file. The 
            files are split into two files to make it easier to read. See 
            output_stats_per_region for details.
    """
    # Uses argparse to organize command line arguments.
    args_to_parse = init_argparse()
    parsed = args_to_parse.parse_args()
    # Iterates through and executes main functions.
    organize_files(parsed.fasta, parsed.kmer,
                    parsed.motif, parsed.title,
                    parsed.output)

# Executes main function.
if __name__ == "__main__":
    main()