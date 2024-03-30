import argparse
from natsort import natsorted

def join_binding_regions(sample_1_bedfile, additional_bedfiles, joined_gtf_file):
    """Takes input BED files and joins all regions into a dictionary,
        removing all duplicate/overlapping regions.

        Arguments:
        sample_1_bedfile -- BED file produced from Piranha or other CLIP
            read pileup tool.
        additional_bedfiles -- List of additional BED files (strings) 
            needed to join first BED file. BED files should be produced 
            from Piranha or other CLIP read pileup tool. Contains one 
            or more bed files.
        joined_gtf_file -- Output GTF file. Contains nonrepetitive,
            comprehensive read pileup regions from each input sample.

        Output:    
        No output. Uses output_file to write out data.
    """
    # Dictionary used to output non-redundant BED entries.
    region_dictionary_out = {}
    # Boolean used to check if first BED file has be iterated through. Used
    # to update region_dictionary with region_dictionary_out for multiple
    # BED file comparisons.
    first_bed_boolean = False
    # Set used to save all output events. Used to output events that did
    # not get joined with a previous event.
    # region_set_added- {"chrom_cord1_cord2"}
    region_set_added = set()
    # Adds all BED lines to a dictionary. 
    # {"chrom_cord1_cord2": "line1_clean"}
    region_dictionary = add_events_to_dictionary(sample_1_bedfile)
    # Iterates through each additional BED file and checks to see if they
    # are already saved in region_dictionary. If so, passes. If one coordinate
    # is present, checks for the larger region and makes sure the larger region
    # is saved in the dictionary.
    for bedfile in additional_bedfiles:
        # Adds all BED lines to a dictionary. This is used to ensure all events
        # get output after comparing all events. This is done by adding joined
        # events to region_set_added.
        # region_dictionary_check- {"chrom_cord1_cord2": "line1_clean"}
        region_dictionary_check = add_events_to_dictionary(bedfile)
        # Updates region_dictionary to with new data if previous BED file added.
        if first_bed_boolean:
            region_dictionary = region_dictionary_out
            # Makes sure the two dictionaries are not pointing to the same object.
            region_dictionary_out = {}
            region_set_added = set()
        for key_2 in region_dictionary_check:
            line_clean = region_dictionary_check.get(key_2)
            line_sep = line_clean.split("_")
            chrom2, cord12, cord22, strand2 = line_sep[0], int(line_sep[1]), int(line_sep[2]), line_sep[5]       
            for key in region_dictionary:
                key_split = key.split("_")
                chrom, cord1, cord2, strand = key_split[0], int(key_split[1]), int(key_split[2]), key_split[3]
                # If already present in dictionary, adds data to output dictionary.
                if ((chrom == chrom2 and strand == strand2) and (cord1 == cord12) and (cord2 == cord22)):
                    region_dictionary_out[key] = region_dictionary.get(key)
                    region_set_added.add(f"{chrom2}_{cord12}_{cord22}_{strand2}")
                    region_set_added.add(f"{chrom}_{cord1}_{cord2}_{strand}")
                # If one of the coordinates matches, checks for larger coordinate. Whatever
                # event is larger gets added to output dictionary.    
                elif (((chrom == chrom2 and strand == strand2) and (cord1 == cord12)) 
                        or ((chrom == chrom2 and strand == strand2) and (cord2 == cord22))):
                    if ((cord1 == cord12) and (cord2 < cord22)):
                        updated_line = f"{chrom}_{cord1}_{cord22}_{line_sep[3]}_{line_sep[4]}_{line_sep[5]}_{line_sep[6]}"
                        region_dictionary_out[f"{chrom}_{cord1}_{cord22}_{strand}"] = updated_line
                    elif ((cord1 == cord12) and (cord2 > cord22)):
                        updated_line = f"{chrom}_{cord1}_{cord2}_{line_sep[3]}_{line_sep[4]}_{line_sep[5]}_{line_sep[6]}"
                        region_dictionary_out[f"{chrom}_{cord1}_{cord2}_{strand}"] = updated_line
                    elif ((cord2 == cord22) and (cord1 < cord12)):
                        updated_line = f"{chrom}_{cord1}_{cord2}_{line_sep[3]}_{line_sep[4]}_{line_sep[5]}_{line_sep[6]}"
                        region_dictionary_out[f"{chrom}_{cord1}_{cord2}_{strand}"] = updated_line
                    elif ((cord2 == cord22) and (cord1 > cord12)):
                        updated_line = f"{chrom}_{cord12}_{cord2}_{line_sep[3]}_{line_sep[4]}_{line_sep[5]}_{line_sep[6]}"
                        region_dictionary_out[f"{chrom}_{cord12}_{cord2}_{strand}"] = updated_line        
                    region_set_added.add(f"{chrom2}_{cord12}_{cord22}_{strand2}")
                    region_set_added.add(f"{chrom}_{cord1}_{cord2}_{strand}")
                # Checks for one of the opposite coordinates are the same. 
                elif (((chrom == chrom2 and strand == strand2)) and ((cord1 == cord22) or (cord2 == cord12))):
                    if (cord1 == cord22):
                        updated_line = f"{chrom2}_{cord12}_{cord2}_{line_sep[3]}_{line_sep[4]}_{line_sep[5]}_{line_sep[6]}"
                        region_dictionary_out[f"{chrom2}_{cord12}_{cord2}_{strand2}"] = updated_line
                    else:
                        updated_line = f"{chrom2}_{cord1}_{cord22}_{line_sep[3]}_{line_sep[4]}_{line_sep[5]}_{line_sep[6]}"
                        region_dictionary_out[f"{chrom2}_{cord1}_{cord22}_{strand2}"] = updated_line
                    region_set_added.add(f"{chrom2}_{cord12}_{cord22}_{strand2}")
                    region_set_added.add(f"{chrom}_{cord1}_{cord2}_{strand}") 
                # Checks if the initial event is overlapping and is larger than the event being compared for both coordinates.
                elif ((chrom == chrom2 and strand == strand2) and (cord1 <= cord12) and (cord2 >= cord22)):
                    updated_line = f"{chrom2}_{cord1}_{cord2}_{line_sep[3]}_{line_sep[4]}_{line_sep[5]}_{line_sep[6]}"
                    region_dictionary_out[f"{chrom2}_{cord1}_{cord2}_{strand2}"] = updated_line
                    region_set_added.add(f"{chrom2}_{cord12}_{cord22}_{strand2}")
                    region_set_added.add(f"{chrom}_{cord1}_{cord2}_{strand}")
                # Checks if the second event is overlapping and is larger than the event being compared both coordinates.
                elif ((chrom == chrom2 and strand == strand2) and (cord1 >= cord12) and (cord2 <= cord22)):
                    updated_line = f"{chrom2}_{cord12}_{cord22}_{line_sep[3]}_{line_sep[4]}_{line_sep[5]}_{line_sep[6]}"
                    region_dictionary_out[f"{chrom2}_{cord12}_{cord22}_{strand2}"] = updated_line
                    region_set_added.add(f"{chrom2}_{cord12}_{cord22}_{strand2}")
                    region_set_added.add(f"{chrom}_{cord1}_{cord2}_{strand}")
                # Checks if events are overlapping but for first event, first coordinate starts earlier and second event, second
                # coordinate ends later in alignment to target sequence.
                elif ((chrom == chrom2 and strand == strand2) and (cord1 <= cord12) and (cord12 <= cord2) and (cord22 >= cord2)):
                    updated_line = f"{chrom2}_{cord1}_{cord22}_{line_sep[3]}_{line_sep[4]}_{line_sep[5]}_{line_sep[6]}"
                    region_dictionary_out[f"{chrom2}_{cord1}_{cord22}_{strand2}"] = updated_line
                    region_set_added.add(f"{chrom2}_{cord12}_{cord22}_{strand2}")
                    region_set_added.add(f"{chrom}_{cord1}_{cord2}_{strand}")
                # Checks if events are overlapping but but second event, first coordinate starts earlier and first event, second
                # coordinate ends later in alignment to target sequence.    
                elif ((chrom == chrom2 and strand == strand2) and (cord1 >= cord12) and (cord1 <= cord22) and (cord2 >= cord22)):
                    updated_line = f"{chrom2}_{cord12}_{cord2}_{line_sep[3]}_{line_sep[4]}_{line_sep[5]}_{line_sep[6]}"
                    region_dictionary_out[f"{chrom2}_{cord12}_{cord2}_{strand2}"] = updated_line
                    region_set_added.add(f"{chrom2}_{cord12}_{cord22}_{strand2}")
                    region_set_added.add(f"{chrom}_{cord1}_{cord2}_{strand}")
        # Sets to true since one comparison has been completed.    
        first_bed_boolean = True
        # Loops through events in original and second input BED file and checks for events that
        # did not get combined and hence not output. Adds these events to dictionary for output.
        for key in region_dictionary_check:
            if key not in region_set_added:
                region_dictionary_out[key] = region_dictionary_check.get(key)
        for key in region_dictionary:
            if key not in region_set_added:
                region_dictionary_out[key] = region_dictionary.get(key)                                                        
    output_file(region_dictionary_out, joined_gtf_file)                            

def last_check_list(sorted_list):
    """Checks list for any missed overlapping events.
    
        Arguments:
        sorted_list -- Naturally sorted list with BED event lines.
    
        Output:    
        region_dictionary -- Filled output dictionary with each read pileup
            region. Contains genomic location information as the value and
            the full bed line as the key. Should not have any overlapping regions.
            {"chrom_cord1_cord2": "line_clean"}
    """
    # Dictionary that stores each read pileup region and outputs
    region_dictionary = {}
    # Saves last event for positive and negative added to dictionary. 
    # Used to check for overlapping events.
    last_event_chrom_plus = ""
    last_event_cord1_plus = ""
    last_event_cord2_plus = ""
    last_event_strand_plus = ""
    last_event_chrom_neg = ""
    last_event_cord1_neg = ""
    last_event_cord2_neg = ""
    # Iterates through each key in dictionary and puts in region_dictionary for
    # output.
    for item in sorted_list:
        line_clean = item
        line_split = line_clean.split("_")
        if line_split[5] == "+":
            # Checks if event is overlapping with previous event. If so, updates the event
            # to include the new coordinate. If not, adds event to dictionary. Also removes
            # last event from dictionary to avoid duplicate events.
            if ((last_event_chrom_plus == line_split[0] and last_event_strand_plus == line_split[5])
                and (last_event_cord2_plus == line_split[1])):
                updated_line = f"{last_event_chrom_plus}_{last_event_cord1_plus}_{line_split[2]}_{line_split[3]}_{line_split[4]}_{line_split[5]}_{line_split[6]}"
                chrom_cord1_cord2_strand = (f"{last_event_chrom_plus}_{last_event_cord1_plus}_{line_split[2]}_{line_split[5]}")
                region_dictionary[chrom_cord1_cord2_strand] = updated_line
                # Removes last event from dictionary. This is done to avoid duplicate events.
                del region_dictionary[f"{last_event_chrom_plus}_{last_event_cord1_plus}_{last_event_cord2_plus}_{last_event_strand_plus}"]
                # Last event cord1 remains the same because it is the start of the event.
                last_event_chrom_plus = last_event_chrom_plus
                last_event_cord1_plus = last_event_cord1_plus
                last_event_cord2_plus = line_split[2]
                last_event_strand_plus = line_split[5]
            # Checks the first two coordinates to see if they are the same, if they are, checks
            # for the second larger coordinate and updates the event to include the larger coordinate.    
            elif ((last_event_chrom_plus == line_split[0] and last_event_strand_plus == line_split[5])
                and (last_event_cord1_plus == line_split[1])):
                if (int(last_event_cord2_plus) > int(line_split[2])):
                    updated_line = f"{last_event_chrom_plus}_{last_event_cord1_plus}_{last_event_cord2_plus}_{last_event_strand_plus}_{line_split[4]}_{line_split[5]}_{line_split[6]}"
                    chrom_cord1_cord2_strand = (f"{last_event_chrom_plus}_{last_event_cord1_plus}_{last_event_cord2_plus}_{last_event_strand_plus}")
                    region_dictionary[chrom_cord1_cord2_strand] = updated_line
                    # Removes last event from dictionary. This is done to avoid duplicate events.
                    del region_dictionary[f"{last_event_chrom_plus}_{last_event_cord1_plus}_{last_event_cord2_plus}_{last_event_strand_plus}"]
                    # Last event cord2 remains the same because it is the end of the event.
                    last_event_chrom_plus = last_event_chrom_plus
                    last_event_cord1_plus = last_event_cord1_plus
                    last_event_cord2_plus = last_event_cord2_plus
                    last_event_strand_plus = last_event_strand_plus
                else:
                    updated_line = f"{last_event_chrom_plus}_{last_event_cord1_plus}_{line_split[2]}_{line_split[3]}_{line_split[4]}_{line_split[5]}_{line_split[6]}"
                    chrom_cord1_cord2_strand = (f"{last_event_chrom_plus}_{last_event_cord1_plus}_{line_split[2]}_{line_split[5]}")
                    region_dictionary[chrom_cord1_cord2_strand] = updated_line
                    # Removes last event from dictionary. This is done to avoid duplicate events.
                    del region_dictionary[f"{last_event_chrom_plus}_{last_event_cord1_plus}_{last_event_cord2_plus}_{last_event_strand_plus}"]
                    # Last event cord2 remains the same because it is the end of the event.
                    last_event_chrom_plus = last_event_chrom_plus
                    last_event_cord1_plus = last_event_cord1_plus
                    last_event_cord2_plus = line_split[2]
                    last_event_strand_plus = line_split[5]
            elif ((last_event_chrom_plus == line_split[0] and last_event_strand_plus == line_split[5])
                and (last_event_cord2_plus == line_split[2])):
                if (int(last_event_cord1_plus) < int(line_split[1])):
                    updated_line = f"{last_event_chrom_plus}_{last_event_cord1_plus}_{last_event_cord2_plus}_{last_event_strand_plus}_{line_split[4]}_{line_split[5]}_{line_split[6]}"
                    chrom_cord1_cord2_strand = (f"{last_event_chrom_plus}_{last_event_cord1_plus}_{last_event_cord2_plus}_{last_event_strand_plus}")
                    region_dictionary[chrom_cord1_cord2_strand] = updated_line
                    # Removes last event from dictionary. This is done to avoid duplicate events.
                    del region_dictionary[f"{last_event_chrom_plus}_{last_event_cord1_plus}_{last_event_cord2_plus}_{last_event_strand_plus}"]
                    # Last event cord1 remains the same because it is the start of the event.
                    last_event_chrom_plus = last_event_chrom_plus
                    last_event_cord1_plus = last_event_cord1_plus
                    last_event_cord2_plus = last_event_cord2_plus
                    last_event_strand_plus = last_event_strand_plus
                else:
                    updated_line = f"{last_event_chrom_plus}_{line_split[1]}_{last_event_cord2_plus}_{line_split[3]}_{line_split[4]}_{line_split[5]}_{line_split[6]}"
                    chrom_cord1_cord2_strand = (f"{last_event_chrom_plus}_{line_split[1]}_{last_event_cord2_plus}_{line_split[5]}")
                    region_dictionary[chrom_cord1_cord2_strand] = updated_line
                    # Removes last event from dictionary. This is done to avoid duplicate events.
                    del region_dictionary[f"{last_event_chrom_plus}_{last_event_cord1_plus}_{last_event_cord2_plus}_{last_event_strand_plus}"]
                    # Last event cord1 remains the same because it is the start of the event.
                    last_event_chrom_plus = last_event_chrom_plus
                    last_event_cord1_plus = line_split[1]
                    last_event_cord2_plus = last_event_cord2_plus
                    last_event_strand_plus = last_event_strand_plus
            else:
                chrom_cord1_cord2_strand = (f"{line_split[0]}_{line_split[1]}_{line_split[2]}_{line_split[5]}")
                region_dictionary[chrom_cord1_cord2_strand] = line_clean
                last_event_chrom_plus = line_split[0]
                last_event_cord1_plus = line_split[1]
                last_event_cord2_plus = line_split[2]
                last_event_strand_plus = line_split[5]
        elif line_split[5] == "-":
            # Checks if event is overlapping with previous event. If so, updates the event
            # to include the new coordinate. If not, adds event to dictionary. Also removes
            # last event from dictionary to avoid duplicate events.
            if ((last_event_chrom_neg == line_split[0] and last_event_strand_neg == line_split[5])
                and (last_event_cord2_neg == line_split[1])):
                updated_line = f"{last_event_chrom_neg}_{last_event_cord1_neg}_{line_split[2]}_{line_split[3]}_{line_split[4]}_{line_split[5]}_{line_split[6]}"
                chrom_cord1_cord2_strand = (f"{last_event_chrom_neg}_{last_event_cord1_neg}_{line_split[2]}_{line_split[5]}")
                region_dictionary[chrom_cord1_cord2_strand] = updated_line
                # Removes last event from dictionary. This is done to avoid duplicate events.
                del region_dictionary[f"{last_event_chrom_neg}_{last_event_cord1_neg}_{last_event_cord2_neg}_{last_event_strand_neg}"]
                # Last event cord1 remains the same because it is the start of the event.
                last_event_chrom_neg = last_event_chrom_neg
                last_event_cord1_neg = last_event_cord1_neg
                last_event_cord2_neg = line_split[2]
                last_event_strand_neg = line_split[5]
            # Checks the first two coordinates to see if they are the same, if they are, checks
            # for the second larger coordinate and updates the event to include the larger coordinate.    
            elif ((last_event_chrom_neg == line_split[0] and last_event_strand_neg == line_split[5])
                and (last_event_cord1_neg == line_split[1])):
                if (int(last_event_cord2_neg) > int(line_split[2])):
                    updated_line = f"{last_event_chrom_neg}_{last_event_cord1_neg}_{last_event_cord2_neg}_{last_event_strand_neg}_{line_split[4]}_{line_split[5]}_{line_split[6]}"
                    chrom_cord1_cord2_strand = (f"{last_event_chrom_neg}_{last_event_cord1_neg}_{last_event_cord2_neg}_{last_event_strand_neg}")
                    region_dictionary[chrom_cord1_cord2_strand] = updated_line
                    # Removes last event from dictionary. This is done to avoid duplicate events.
                    del region_dictionary[f"{last_event_chrom_neg}_{last_event_cord1_neg}_{last_event_cord2_neg}_{last_event_strand_neg}"]
                    # Last event cord2 remains the same because it is the end of the event.
                    last_event_chrom_neg = last_event_chrom_neg
                    last_event_cord1_neg = last_event_cord1_neg
                    last_event_cord2_neg = last_event_cord2_neg
                    last_event_strand_neg = last_event_strand_neg
                else:
                    updated_line = f"{last_event_chrom_neg}_{last_event_cord1_neg}_{line_split[2]}_{line_split[3]}_{line_split[4]}_{line_split[5]}_{line_split[6]}"
                    chrom_cord1_cord2_strand = (f"{last_event_chrom_neg}_{last_event_cord1_neg}_{line_split[2]}_{line_split[5]}")
                    region_dictionary[chrom_cord1_cord2_strand] = updated_line
                    # Removes last event from dictionary. This is done to avoid duplicate events.
                    del region_dictionary[f"{last_event_chrom_neg}_{last_event_cord1_neg}_{last_event_cord2_neg}_{last_event_strand_neg}"]
                    # Last event cord2 remains the same because it is the end of the event.
                    last_event_chrom_neg = last_event_chrom_neg
                    last_event_cord1_neg = last_event_cord1_neg
                    last_event_cord2_neg = line_split[2]
                    last_event_strand_neg = line_split[5]
            elif ((last_event_chrom_neg == line_split[0] and last_event_strand_neg == line_split[5])
                and (last_event_cord2_neg == line_split[2])):
                if (int(last_event_cord1_neg) < int(line_split[1])):
                    updated_line = f"{last_event_chrom_neg}_{last_event_cord1_neg}_{last_event_cord2_neg}_{last_event_strand_neg}_{line_split[4]}_{line_split[5]}_{line_split[6]}"
                    chrom_cord1_cord2_strand = (f"{last_event_chrom_neg}_{last_event_cord1_neg}_{last_event_cord2_neg}_{last_event_strand_neg}")
                    region_dictionary[chrom_cord1_cord2_strand] = updated_line
                    # Removes last event from dictionary. This is done to avoid duplicate events.
                    del region_dictionary[f"{last_event_chrom_neg}_{last_event_cord1_neg}_{last_event_cord2_neg}_{last_event_strand_neg}"]
                    # Last event cord1 remains the same because it is the start of the event.
                    last_event_chrom_neg = last_event_chrom_neg
                    last_event_cord1_neg = last_event_cord1_neg
                    last_event_cord2_neg = last_event_cord2_neg
                    last_event_strand_neg = last_event_strand_neg
                else:
                    updated_line = f"{last_event_chrom_neg}_{line_split[1]}_{last_event_cord2_neg}_{line_split[3]}_{line_split[4]}_{line_split[5]}_{line_split[6]}"
                    chrom_cord1_cord2_strand = (f"{last_event_chrom_neg}_{line_split[1]}_{last_event_cord2_neg}_{line_split[5]}")
                    region_dictionary[chrom_cord1_cord2_strand] = updated_line
                    # Removes last event from dictionary. This is done to avoid duplicate events.
                    del region_dictionary[f"{last_event_chrom_neg}_{last_event_cord1_neg}_{last_event_cord2_neg}_{last_event_strand_neg}"]
                    # Last event cord1 remains the same because it is the start of the event.
                    last_event_chrom_neg = last_event_chrom_neg
                    last_event_cord1_neg = line_split[1]
                    last_event_cord2_neg = last_event_cord2_neg
                    last_event_strand_neg = last_event_strand_neg
            else:
                chrom_cord1_cord2_strand = (f"{line_split[0]}_{line_split[1]}_{line_split[2]}_{line_split[5]}")
                region_dictionary[chrom_cord1_cord2_strand] = line_clean
                last_event_chrom_neg = line_split[0]
                last_event_cord1_neg = line_split[1]
                last_event_cord2_neg = line_split[2]
                last_event_strand_neg = line_split[5]                                                                                                                         
    return region_dictionary

def add_events_to_dictionary(sample_1_bedfile):
    """Takes an input BED file and adds all regions into a dictionary.

        Arguments:
        sample_1_bedfile -- BED file produced from Piranha or other CLIP
            read pileup tool.

        Output:    
        region_dictionary -- Filled output dictionary with each read pileup
            region. Contains genomic location information as the value and
            the full bed line as the key. {"chrom_cord1_cord2_strand": "line_clean"}.
    """                      
    # Dictionary that stores each read pileup region.
    region_dictionary = {}
    # Saves last event added to dictionary. Used to check for overlapping events.
    last_event_chrom = ""
    last_event_cord1 = ""
    last_event_cord2 = ""
    last_event_strand = ""
    # Iterates through each line in the BED file and adds it to the dictionary.                         
    with open(sample_1_bedfile) as open_bedfile:
        for line in open_bedfile:
            # Replaces tabs with underscore. This makes it so the literal tab
            # characters do not get stored in the dictionary, making it more difficult
            # for a clean output.
            line_clean = line.strip("\n").replace("\t", "_")
            line_split = line_clean.split("_")
            # Checks if event is overlapping with previous event. If so, updates the event
            # to include the new coordinate. If not, adds event to dictionary. Also removes
            # last event from dictionary to avoid duplicate events.
            if ((last_event_chrom == line_split[0] and last_event_strand == line_split[5])
                and (last_event_cord2 == line_split[1])):
                updated_line = f"{last_event_chrom}_{last_event_cord1}_{line_split[2]}_{line_split[3]}_{line_split[4]}_{line_split[5]}_{line_split[6]}"
                chrom_cord1_cord2_strand = (f"{last_event_chrom}_{last_event_cord1}_{line_split[2]}_{line_split[5]}")
                region_dictionary[chrom_cord1_cord2_strand] = updated_line
                # Removes last event from dictionary. This is done to avoid duplicate events.
                del region_dictionary[f"{last_event_chrom}_{last_event_cord1}_{last_event_cord2}_{last_event_strand}"]
                # Last event cord1 remains the same because it is the start of the event.
                last_event_chrom = last_event_chrom
                last_event_cord1 = last_event_cord1
                last_event_cord2 = line_split[2]
                last_event_strand = line_split[5]
            else:
                chrom_cord1_cord2_strand = (f"{line_split[0]}_{line_split[1]}_{line_split[2]}_{line_split[5]}")
                region_dictionary[chrom_cord1_cord2_strand] = line_clean
                last_event_chrom = line_split[0]
                last_event_cord1 = line_split[1]
                last_event_cord2 = line_split[2]
                last_event_strand = line_split[5]
    return region_dictionary

def output_file(region_dictionary_out, joined_gtf_file):
    """Takes dictionary of unique read pileup regions and writes each region
        out into a TSV file.

        Arguments:
        region_dictionary_out -- Filled output dictionary with each read pileup
            region. Contains genomic location information as the value and
            the full bed line as the key. {"chrom_cord1_cord2": "line_clean"}

        Output:    
        joined_gtf_file -- Output GTF file. Contains nonrepetitive, comprehensive
            read pileup regions from each input sample
    """
    # Opens joined_bed_file to write out data, adds value to list, sorts list,
    # and writes out each entry.   
    joined_gtf_file_out = open(joined_gtf_file, 'w')
    # Initializes list used for sorting events before output.
    sorted_list = []
    # Iterates through dictionary and adds event to list.                       
    for key in region_dictionary_out:
        sorted_list.append(region_dictionary_out.get(key))
    # Sorts list.    
    sorted_list = natsorted(sorted_list)
    # Performs last check to make sure there are no 
    #overlapping events in the output dictionary.
    region_dictionary = last_check_list(sorted_list)
    # Makes empty list to populate from scratch.
    sorted_list = []
    # Iterates through dictionary and adds deduplicated events to list again
    # for output.                       
    for key in region_dictionary:
        sorted_list.append(region_dictionary_out.get(key)) 
    # Outputs sorted list in GTF format. Outputs gene and transcript annotation for each event.
    for item in sorted_list:
        item_split = item.split("_")
        # Changes first coordinate to 1 if 0. This is done to avoid errors in downstream
        # analysis.
        if item_split[1] == "0":
            coordinate_1 = "1"
        else:
            coordinate_1 = item_split[1]    
        joined_gtf_file_out.write(item_split[0] + "\t" + "pir" + "\t" + "gene" + "\t" + coordinate_1 + "\t" + item_split[2] \
                                  + "\t" + item_split[4] + "\t" + item_split[5] + "\t" + item_split[6] + "\t" + 'gene_id "' \
                                  + item_split[0] + "~" + coordinate_1 + "~" + item_split[2] + "~" + item_split[5] +'g"; ' \
                                  + 'gene_version "1"; ' + "\n")
        joined_gtf_file_out.write(item_split[0] + "\t" + "pir" +"\t"+ "transcript" + "\t" + coordinate_1 + "\t" + item_split[2] \
                                  + "\t" + item_split[4] + "\t" + item_split[5] + "\t" + item_split[6] + "\t" + 'gene_id "' \
                                  + item_split[0] + "~" + coordinate_1 + "~" + item_split[2] + "~"  +item_split[5] + '"; ' \
                                  + 'gene_version "1"; ' + 'transcript_id "' +item_split[0] + "~" + coordinate_1 + "~" \
                                  + item_split[2] + "~" + item_split[5] + 't"; ' + "\n")
    joined_gtf_file_out.close()

def init_argparse():
    """Initiates the use of argparse. Returns parsed
        command line arguments.

        Arguments:
        See command line argument descriptions below.

        Output:
        Returns parser object. Object is used for function input.
    """
    # Script description.
    parser = argparse.ArgumentParser(description = "Joins regions of read pileups produced from \
                                                    Piranha or other CLIP read pileup tool.")        
    # Command line arguments and descriptions.
    parser.add_argument("-b1", "--sample_1_bedfile", action = "store", type = str, 
                        help="BED file produced from Piranha or other CLIP read pileup tool.",
                        required = True)
    parser.add_argument("-b2", "--additional_bedfiles", action = "store", type = str, nargs='+', 
                        help = "Additional BED files needed to join first BED file. BED files \
                        should be produced from Piranha or other CLIP read pileup tool. Add \
                        each BED file needed to be joined separated by a space.", 
                        required = True)
    parser.add_argument("-o", "--joined_gtf_file", action = "store", type = str, 
                        help="Output GTF file.", 
                        required = False)                                                                              
    return parser

def main():
    """Main script function. Executes other functions.

        Arguments:
        sample_1_bedfile -- BED file produced from Piranha or other CLIP
            read pileup tool.
        additional_bedfiles -- Additional BED files needed to join first BED file.
            BED files should be produced from Piranha or other CLIP read pileup tool.
            This is a list of one or more BED files.

        Output:    
        joined_gtf_file -- Output GTF file. Contains nonrepetitive, comprehensive
            read pileup regions from each input sample.
    """
    # Uses argparse to organize command line arguments.
    args_to_parse = init_argparse()
    parsed = args_to_parse.parse_args()
    # Iterates through and executes main functions.
    join_binding_regions(parsed.sample_1_bedfile, parsed.additional_bedfiles,
                         parsed.joined_gtf_file)

# Executes main function.
if __name__ == "__main__":
    main()