import argparse

def make_splicing_dictionary(rmats_file, output):
    """Generates a dictionary of splicing events.

        Arguments:
        rmats_file -- rMATS alternative skipped exon splicing 
            output file (SE.MATS.JCEC.txt file). Should have 
            been generated using the same target organism 
            GTF file. All other types of alternative splicing 
            can be input except for retained intron (different 
            number of columns).
        output -- String prefix for output file. Should contain
            the full path to the output directory if a different
            directory is being used.

        Output:
        splicing_dictionary -- Dictionary of splicing events
            with exon coordinates as keys and additional
            splicing information as the value.
            {"gene-chromosome-exon_start-exon_end-strand":
             [ensemble_id, up-start, up-end, down-start, down-end, 
              delta_psi, regulation]}
        filtered_splicing_out -- rMATS style output file that is
            filtered by delta PSI and an FDR of 10% or less.
            File name- output + "filtered_SE.MATS.JCEC.txt".  
    """
    # Loops through splicing file and grabs exon coordinates.      
    rmats_file_open = open(rmats_file, 'r')
    rmats_output = open(output + "filtered_SE.MATS.JCEC.txt", 'w')
    splicing_dictionary = {}
    for line in rmats_file_open:
        # Skips header line. Write out header line to filtered file.
        if line[0] != "I":
            # Prepares variables for dictionary.
            line_pre = line.replace('"' , '').strip("\n")
            line_split = line_pre.split("\t")
            # Filters for a delta PSI of greater than 10% and FDR of less than 10%.
            if (abs(float(line_split[22])) > 0.10 and (float(line_split[19]) < 0.10)):
                # exon_key = [gene, chromosome, exon-start, exon-end, strand]
                exon_key = f"{line_split[2]}~{line_split[3]}~{line_split[5]}~{line_split[6]}~{line_split[4]}"
                ensemble_id = line_split[1].replace('"' , '')
                # Indicates if exon is included or excluded ("Inclusion" or "Exclusion").
                # NOTE- This may vary per experiment and controls. Please analyze your
                # splicing file to determine what positive and negative values indicate.
                regulation = ""
                if float(line_split[22]) > 0: 
                    regulation = "Exclusion"
                else:
                    regulation = "Inclusion"    
                # exon_list = [ensemble_id, up-start, up-end, down-start, down-end, 
                #               delta_psi, regulation]
                exon_list = [ensemble_id, line_split[7], line_split[8], line_split[9], 
                              line_split[10], line_split[22], regulation]
                splicing_dictionary[exon_key] = exon_list
                # Writes filtered splicing events to file.
                rmats_output.write(line)
        else:
            rmats_output.write(line)
    return splicing_dictionary

def match_binding_to_splicing(binding_regions, splicing_dictionary, distance, 
                              output, region_csv):
    """Matches binding regions to splicing events.
    
        Arguments:
        binding_regions -- DESeq2 output file with binding
            regions. Should be annotated and filtered file,
            although any DESeq2 output file can be used.
        splicing_dictionary -- Dictionary of splicing events
            with exon coordinates as keys and additional
            splicing information as the value.
            {"gene-chromosome-exon_start-exon_end-strand":
             [ensemble_id, up-start, up-end, down-start, down-end, 
              delta_psi, regulation]}
        distance -- Integer distance skipped exon outer
            cordinates should be to be considered within
            regulated distance of binding region.
        output -- Prefix for output file.
        region_csv -- CSV file of random motif enrichment values
            for each binding region from all input files. Was generated
            from the random_enrichment_motify.py script. File name-
            output_file_prefix + "_motifs_per_region.csv". Format -
            Region Name, Random Motif Enrichment Value, Normalized.  

        Output:
        splicing_binding_out -- CSV file with splicing and binding
            information for each region.
    """
    # Opens binding regions file and creates output file.
    binding_regions_open = open(binding_regions, 'r')
    splicing_binding_out = open(output + "splicing_binding.csv", 'w')
    # Writes header line to output file.
    splicing_binding_out.write("Region_Name,Gene,Chromosome,Exon-Start,Exon-End,Strand,Ensemble-ID,")
    splicing_binding_out.write("Delta_PSI,Regulation_Type,Binding_Region_Start,Binding_Region_End,")
    splicing_binding_out.write("Normalized_Motif_Enrichment,Region Size,")
    splicing_binding_out.write("Region_Position-(+up)-(-down)-(0in),Base-Mean,Canonical_Splicing_Regulation,")
    splicing_binding_out.write("Normalized_Count\n")
    # Puts regions motif values into a dictionary.
    region_dictionary = collect_regions(region_csv)
    # Loops through binding regions file and matches splicing events.
    for line in binding_regions_open:
        line_split = line.split(",")
        # Grabs normalized counts from the last column for each region.
        normalized_count = line_split[-1].strip("\n")
        # Skips header line for DESeq2 output file.
        if (line.find("~") != -1):
            region_pre = line_split[0].replace('"','').split("~")
            region_size = int(region_pre[2]) - int(region_pre[1])
            # region_key = [chromosome, region-start, region-end, strand, region-size, base-mean]
            region_list = [region_pre[0], region_pre[1], region_pre[2], region_pre[3][0], region_size, line_split[1]]
            for key_pre in splicing_dictionary:
                key = key_pre.split("~")
                # Checks if chromosome and strand match.
                if (key[1] == region_list[0] and key[4] == region_list[3]):
                    # Checks if binding region is within distance of splicing event.
                    if ((abs(int(key[2]) - int(region_list[1])) <= distance) or
                        (abs(int(key[3]) - int(region_list[1])) <= distance) or
                        (abs(int(key[2]) - int(region_list[2])) <= distance) or
                        (abs(int(key[3]) - int(region_list[2])) <= distance)):
                        # Region position is the distance from the start of the region to the
                        # start of the splicing event. A positive region position means the
                        # binding region is 3' (up in coordinate number), a negative number
                        # indicates 5', and 0 indicates inside the exon. Depending on where 
                        # this position is, the canonical splicing regulation will be determined. 
                        # NOTE- This may vary per experiment, RBP, and controls. This is default 
                        # for MBNL1. If region is 5' or in the exon and the splicing event is 
                        # being excluded, the splicing event is canonical. If the region is 3' 
                        # and the splicing event is being included, the splicing event is canonical.
                        region_position = 0
                        canonical_splicing = ""
                        splicing_value = splicing_dictionary[key_pre]
                        # This indicates the region is 3'.
                        if ((int(key[3]) - int(region_list[1]) < 0)):
                            region_position = int(key[3]) - int(region_list[1])
                            if (splicing_value[6] == "Inclusion"):
                                canonical_splicing = "Canonical"
                            else:
                                canonical_splicing = "Non-Canonical"
                        # This indicates the region is 5'.        
                        elif ((int(key[2]) - int(region_list[2]) > 0)):
                            region_position = int(key[2]) - int(region_list[2])
                            if (splicing_value[6] == "Exclusion"):
                                canonical_splicing = "Canonical"
                            else:
                                canonical_splicing = "Non-Canonical"
                        # This indicates the region is in the exon.
                        elif (((int(key[2]) - int(region_list[1]) <= 0) and
                              (int(key[3]) - int(region_list[1]) >= 0)) or
                              ((int(key[2]) - int(region_list[2]) <= 0) and
                                (int(key[3]) - int(region_list[2]) >= 0))):
                            region_position = 0
                            if (splicing_value[6] == "Exclusion"):
                                canonical_splicing = "Canonical"
                            else:
                                canonical_splicing = "Non-Canonical"
                        # Grabs motif values for region.
                        motif_value = 0.0
                        motif_values_list = region_dictionary.get(line_split[0])
                        if motif_values_list != None:
                            motif_value = motif_values_list[0]
                        else:
                            motif_value = 0.0
                        # Writes out splicing and binding information to output file.
                        splicing_binding_out.write(line_split[0] + "," + key[0] + "," + key[1] + "," + 
                                                   key[2] + "," + key[3] + "," + key[4] + "," + splicing_value[0] + ",")
                        splicing_binding_out.write(splicing_value[5] + "," + splicing_value[6] + "," + 
                                                   region_list[1] + "," + region_list[2] + ",")
                        splicing_binding_out.write(str(motif_value) + "," + str(region_list[4])+ ",")
                        print(motif_value)
                        splicing_binding_out.write(str(region_position) + "," + str(region_list[5]) + "," + 
                                                   canonical_splicing + "," + normalized_count + "\n")

def collect_regions(region_csv):
    """Collects binding region motif values into a dictionary.

        Arguments:
        region_csv -- CSV file of random motif enrichment values
            for each binding region from all input files. Was generated
            from the random_enrichment_motify.py script. File name-
            output_file_prefix + "_motifs_per_region.csv". Format -
            Region Name, Random Motif Enrichment Value, Normalized.  

        Output:
        region_dictionary -- Dictionary of binding regions with
            region name as key and motif values as value.
            {region_name: [motif enrichment, normalized]}
    """
    # Opens region CSV file and creates output dictionary.
    region_csv_open = open(region_csv, 'r')
    region_dictionary = {}
    # Skips header line.
    for line in region_csv_open:
        if line[0] != "R":
            line_split = line.strip("\n").split(",")
            # Formats region name to match DeSeq2 output.
            region_pre = line_split[0].split("(")
            strand = region_pre[1].replace(")", "")
            region_name = region_pre[0].replace("-", "~")
            region_name = region_name.replace(":", "~")
            region_name = '"' + region_name + "~" + strand + "g" + '"'
            region_dictionary[region_name] = [line_split[1], line_split[2]]                                         
    return region_dictionary

def init_argparse():
    """Initiates the use of argparse. Returns parsed
        command line arguments.

        Arguments:
        See command line argument descriptions below.

        Output:
        Returns parser object. Object is used for function input.
    """
    # Script description.
    parser = argparse.ArgumentParser(description = "Matches alternative skipped exon splicing events \
                                                    to binding regions within a specified distance.")
    # Command line arguments and descriptions.
    parser.add_argument("-r", "--rmats_file", action = "store", type = str,
                        help="rMATS alternative skipped exon splicing output file \
                              (SE.MATS.JCEC.txt file). Should have been generated \
                              using the same target organism GTF file.", 
                        required = True)
    parser.add_argument("-b", "--binding_regions", action = "store", type = str,
                        help="DESeq2 output file with binding regions. Should be \
                              annotated and filtered file.",
                        required = True)
    parser.add_argument("-m", "--region_csv", action = "store", type = str,
                        help="CSV file of random motif enrichment values for each \
                             binding region from all input files.",
                        required = True)
    parser.add_argument("-d", "--distance", action = "store", type = int,
                        help="Distance skipped exon outer cordinates should be to be \
                            considered within regulated distance of binding region.",
                        required = True)
    parser.add_argument("-o", "--output", action = "store", type = str,
                        help="Prefix for output files.",
                        required = True)                                                                 
    return parser

def main():
    """Main script function. Executes other functions.

        Arguments:
        rmat_file -- rMATS alternative skipped exon splicing 
            output file (SE.MATS.JCEC.txt file). Should have 
            been generated using the same target organism 
            GTF file. All other types of alternative splicing 
            can be input except for retained intron (different 
            number of columns).
        binding_regions -- DESeq2 output file with binding
            regions. Should be annotated and filtered file,
            although any DESeq2 output file can be used.
        region_csv -- CSV file of random motif enrichment values
            for each binding region from all input files. Was generated
            from the random_enrichment_motify.py script. File name-
            output_file_prefix + "_motifs_per_region.csv". Format -
            Region Name, Random Motif Enrichment Value, Normalized.    
        distance -- Integer distance skipped exon outer
            cordinates should be to be considered within
            regulated distance of binding region.    
        output -- Prefix for output file.                   

        Output: 
        filtered_splicing_out -- rMATS style output file that is
            filtered by delta PSI and an FDR of 10% or less.
            File name- output + "filtered_SE.MATS.JCEC.txt".    
        splicing_binding_out -- CSV file with splicing and binding
            information for each region.

    """
    # Uses argparse to organize command line arguments.
    args_to_parse = init_argparse()
    parsed = args_to_parse.parse_args()
    # Iterates through and executes main functions.
    splicing_dictionary = make_splicing_dictionary(parsed.rmats_file, parsed.output)
    match_binding_to_splicing(parsed.binding_regions, splicing_dictionary, parsed.distance, 
                              parsed.output, parsed.region_csv)

# Executes main function.
if __name__ == "__main__":
    main()