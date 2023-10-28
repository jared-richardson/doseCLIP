import argparse

def motif_detector(region_fasta_list, kmer_size, motif_list = ["YGCY"], 
                   motif_title = "YGCY"):
    """Analyzes a region FASTA file for the presence of a motif.

        Arguments:
        region_fasta_list -- List of FASTA files of RBP binding regions.
            Should include previously generated random regions from
            make_random_bed_for_motif.py. The random files are detected
            by searching for "random" in the file name.
        kmer_size -- Size of kmer to be analyzed.
        motif_list -- List of motifs to be analyzed. Default is YGCY.
        motif_title -- Title for motif output file. Default is YGCY.

        Output:    
        file_dictionary -- Dictionary of motif counts for each file.
            {"file_name": {"motif": normalized_count(float)}}
        region_dictionary -- Dictionary of motif counts for each region.
            {"event_name": normalized_motif_count(float)}
        motif_list -- List of motifs to be analyzed. Default is YGCY.    
        motif_title -- Title for motif output file. Default is YGCY.    
    """
    # Used to stores YGCY motif counts for each region for all files.
    # {"event_name": normalized_motif_count(float)}
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
    # Used to store motif counts for each file.
    # {"file_name": {"motif": normalized_count(float)}}
    file_dictionary = {}
    # Used to store name of binding region.
    event_name = ""
    # Loops through each FASTA file and generates normalized
    # motif counts for each file and for each region.
    for file in region_fasta_list:
        # Used to store motif counts.
        # {"motif": normalized_count(float)}
        motif_dictionary = {}
        # Counts total nucleotides. Used for normalization.
        nucleotide_count = 0
        open_file = open(file, 'r')
        # Looks for "/" in directory structure. If found, splits
        # and grabs last element. Otherwise, just uses file name.
        if file.find("/") != -1:
            file_split = file.split("/")
            file_name = file_split[-1]
        else:
            file_name = file 
        # Loops through each line of the FASTA file. If title line
        # title is saved, if nucleotides, motif counts are generated.        
        for line in open_file:
            # Initializes region_motif_dictionary for each region.
            region_motif_dictionary = {}
            # Used to store nucleotide counts for each region.
            nucleotide_count_region = 0
            line_clean = line.strip("\n")
            # Used to detect title line.
            if line[0] != ">":
                # Loops through each nucleotide and generates motif counts. 
                # Generates normalized counts per region and per file.
                for start in range(len(line_clean) - kmer_size + 1):
                    kmer = line_clean[start:start + kmer_size]
                    current_count = motif_dictionary.get(kmer, 0)
                    current_count_region = region_motif_dictionary.get(kmer, 0)
                    motif_dictionary[kmer] = current_count + 1
                    region_motif_dictionary[kmer] = current_count_region + 1 
                nucleotide_count = (nucleotide_count + len(line_clean))
                nucleotide_count_region = (nucleotide_count_region + len(line_clean))
                # Looks for YGCY variants and adds them together, normalizes for 
                # nucleotides, and adds to region_dictionary. If separate motif
                # list is input, uses that instead. Skips random regions from 
                # random region files.
                if file_name.find("random") == -1:
                    # Used to total motif counts for each region.
                    motif_total = 0
                    for kmer in region_motif_dictionary:
                        if (motif_title == "YGCY") and (kmer in motif_list):
                            motif_total = (motif_total + region_motif_dictionary[kmer])
                        elif kmer in motif_list:
                            motif_total = (motif_total + region_motif_dictionary[kmer])
                    if event_name not in region_dictionary:
                        region_dictionary[event_name] = (motif_total / float(nucleotide_count_region))
            else:
                event_name = line_clean.replace(">", "")
        # Loops through motif_dictionary and normalizes counts using nucleotide_count from
        # entire file. Adds to file_dictionary.
        for kmer in motif_dictionary:
            normalized_count = (motif_dictionary[kmer] / float(nucleotide_count))
            motif_dictionary[kmer] = normalized_count        
        file_dictionary[file_name] = motif_dictionary            
        open_file.close() 
    return file_dictionary, region_dictionary, motif_list, motif_title       
            
def output_motifs(file_dictionary, region_dictionary, output_file_prefix,
                  motif_list, motif_title):
    """Analyzes a region FASTA file for the presence of a motif.

        Arguments:
        file_dictionary -- Dictionary of motif counts for each file.
            {"file_name": {"motif": normalized_count(float)}}
        region_dictionary -- Dictionary of motif counts for each region.
            {"event_name": normalized_motif_count(float)}
        output_file_prefix -- Prefix for output files.
        motif_list -- List of motifs to be analyzed. Default is YGCY.
        motif_title -- Title for motif output file. Default is YGCY.

        Output:    
        Output region CSV -- CSV files of random motif enrichment values
            for each binding region from all input files. File name-
            output_file_prefix + "_motifs_per_region.csv". Format -
            Region Name, Motif Per 100 Nucleotides Value. 
        Output motif CSV -- CSV file of random motif enrichment values
            per target motif per file. Last section indicates whether
            motif was able to be normalized using the random files.
            All should be normalized but was added for edge cases.
            File name- output_file_prefix + "total_motifs.csv". 
            Format -
            File Name, Motif, Random Motif Enrichment Value, Normalized.
        Output all motifs CSV -- CSV file of random motif enrichment values
            for all motifs per file. Last section indicates whether
            motif was able to be normalized using the random files.
            All should be normalized but was added for edge cases.
            File name- output_file_prefix + "all_motifs.csv". 
            Format -
            File Name, Motif, Random Motif Enrichment Value, Normalized.    

    """
    # Names and opens output files for writing.
    output_region = (output_file_prefix + "_motifs_per_region.csv")
    output_region_out = open(output_region, 'w')
    output_motif = (output_file_prefix + "_total_motifs.csv")
    output_motif_out = open(output_motif, 'w')
    output_all_motifs = (output_file_prefix + "_all_motifs.csv")
    output_all_motifs_out = open(output_all_motifs, 'w')
    # Writes headers for output files.
    output_region_out.write(f"Region Name, {motif_title} Motif Per 100 Nucleotides Value\n")
    output_motif_out.write(f"File Name, {motif_title} Motif, Random Motif Enrichment Value, Normalized\n")
    output_all_motifs_out.write("File Name, Motif, Random Motif Enrichment Value, Normalized\n")
    # Loops through region_dictionary and outputs region name and
    # motif count per 100 nucleotides (motif per nucleotides is 
    # multiplied by 100).
    for region in region_dictionary:
        output_region_out.write(region + "," + 
                                str(region_dictionary[region] * 100) + "\n")
    random_dictionary = get_average_random(file_dictionary)    
    # Loops through motif_dictionary and outputs file name, motif, and
    # normalized motif count. Does not use random files.
    for file in file_dictionary:
        # String used to indicate normalization.
        normalized = ""
        if file.find("random") == -1:
            for motif in file_dictionary[file]:
                # Uses random_dictionary to normalize motif counts. Checks to see if motif key
                # is in random_dictionary. If not, skips normalization and indicates with string.
                if (random_dictionary.get(motif) != None):
                    motif_normalized = (file_dictionary[file][motif] / random_dictionary[motif])
                    normalized = "Yes"
                else:
                    motif_normalized = (file_dictionary[file][motif])
                    normalized = "No"   
                output_all_motifs_out.write(file + "," + motif + "," + str(motif_normalized) 
                                            + "," + normalized + "\n")
                # Filters using motif_list to output only motifs of interest.
                if motif in motif_list:
                    output_motif_out.write(file + "," + motif + "," + str(motif_normalized)
                                            + "," + normalized + "\n")

def get_average_random(file_dictionary):
    """
    Arguments:
    file_dictionary -- Dictionary of motif counts for each file.
        {"file_name": {"motif": normalized_count(float)}}

    Output:
    random_dictionary -- Dictionary of average motif counts for
        all random files.
        {"motif": normalized_count(float)}
    """
    # Used to store average random motif counts.
    random_dictionary = {}
    # Used to determine total random file count.
    random_file_count = 0
    # Loops through motifs, pulls out random files, and sums of all
    # random motif counts.
    for file in file_dictionary:
        if file.find("random") != -1:
            random_file_count += 1
            for motif in file_dictionary[file]:
                current_count = random_dictionary.get(motif, 0)
                random_dictionary[motif] = current_count + file_dictionary[file][motif]
    # Loops through random_dictionary and divides by total number of 
    # random files to get average, if more than one random file.
    if random_file_count > 1:           
        for motif in random_dictionary:
            random_dictionary[motif] = (random_dictionary[motif] / float(random_file_count))
    return random_dictionary                   

def init_argparse():
    """Initiates the use of argparse. Returns parsed
        command line arguments.

        Arguments:
        See command line argument descriptions below.

        Output:
        Returns parser object. Object is used for function input.
    """
    # Script description.
    parser = argparse.ArgumentParser(description = "Gets normalized random motif enrichment values for each \
                                     file input and each region in input files.")        
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
        Output region CSV -- CSV files of random motif enrichment values
            for each binding region from all input files. File name-
            output_file_prefix + "_motifs_per_region.csv". Format -
            Region Name, Motif Per 100 Nucleotides Value. 
        Output motif CSV -- CSV file of random motif enrichment values
            per target motif per file. Last section indicates whether
            motif was able to be normalized using the random files.
            All should be normalized but was added for edge cases.
            File name- output_file_prefix + "total_motifs.csv". 
            Format -
            File Name, Motif, Random Motif Enrichment Value, Normalized.
        Output all motifs CSV -- CSV file of random motif enrichment values
            for all motifs per file. Last section indicates whether
            motif was able to be normalized using the random files.
            All should be normalized but was added for edge cases.
            File name- output_file_prefix + "all_motifs.csv". 
            Format -
            File Name, Motif, Random Motif Enrichment Value, Normalized. 

    """
    # Uses argparse to organize command line arguments.
    args_to_parse = init_argparse()
    parsed = args_to_parse.parse_args()
    # Iterates through and executes main functions.
    (file_dictionary, region_dictionary, 
     motif_list, motif_title) = motif_detector(parsed.fasta, parsed.kmer,
                                              parsed.motif, parsed.title)
    output_motifs(file_dictionary, region_dictionary, parsed.output, 
                  motif_list, motif_title) 

# Executes main function.
if __name__ == "__main__":
    main()