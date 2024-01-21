import argparse

def get_motif_enrichment(fastq_list, fastq_sm_list, kmer_size):
    """Analyzes sets of paired-end FASTQ files with SM Inputs 
        for motif enrichment, RBNS style.

        Arguments:
        fastq_list -- (-f) List of paired-end FASTQ files of doseCLIP
            samples. Should contain all replicates and read1 and 
            read2 should be in same order.
        fastq_sm_list -- (-fs) List of paired-end FASTQ files of SM samples 
            pertaining to each and all input doseCLIP files. Should
            be in same order as doseCLIP samples, including with read1
            and read2. 
        kmer_size -- (-k) Size of kmer to be analyzed.

        Output:    
        final_dictionary -- Dictionary of average motif enrichment
            values using all input files. 
            {"motif": normalized_count(float)}.
        file_dictionary -- Dictionary of motif counts for each file.
            {"file_name": {"motif": normalized_count(float)}}.    
    """
    # Used to store motif counts for each file.
    # {"file_name": {"motif": normalized_count(float)}}
    file_dictionary = {}
    # Used to store final motif enrichment values from all files.
    # {"motif": normalized_motif_count(float)}
    final_dictionary = {}
    # Loops through each FASTQ file and generates normalized
    # motif counts for each file and for each region.
    for file in fastq_list:
        # Used to grab SM input file.
        index = fastq_list.index(file)
        file_sm = fastq_sm_list[index]
        # Gets file names and motif dictionaries.
        file_name = get_file_name(file)
        file_name_sm = get_file_name(file_sm)
        motif_dictionary = motif_detector(file, kmer_size)
        motif_dictionary_sm = motif_detector(file_sm, kmer_size)
        # Normalizes motif counts for each file using SM input.
        normalized_dictionary = normalize_motifs(motif_dictionary, motif_dictionary_sm)
        file_dictionary[file_name + "-" + file_name_sm] = normalized_dictionary
    # Dictionary used to count if a motif was present in each files.
    # Used for averaging enrichment values.
    count_dictionary = {}        
    # Loops through file_dictionary and generates average motif enrichment
    # values using all files.
    for file_name_full in file_dictionary:
        # Loops through each motif in each file.
        for motif in file_dictionary[file_name_full]:
            # Used to store average motif enrichment values.
            current_value = final_dictionary.get(motif, 0)
            final_dictionary[motif] = (current_value + 
                                       file_dictionary[file_name_full][motif])
            # Used to count if motif was present in each file.
            current_count = count_dictionary.get(motif, 0)
            count_dictionary[motif] = (current_count + 1)
    # Uses count_dictionary to average motif enrichment values.
    for motif in final_dictionary:
        final_dictionary[motif] = (final_dictionary[motif] / 
                                   float(count_dictionary[motif]))  
    return (final_dictionary, file_dictionary)   
            
def motif_detector(fastq_file, kmer_size):
    """Gets motif enrichment values per FASTQ file.

    Arguments:
    fastq_file -- FASTQ file to be analyzed for motifs.
    kmer_size -- Size of kmer to be analyzed.

    Output:
    motif_dictionary -- Dictionary of motif counts for each file.
        {"motif": normalized_count(float)}
    """
    # Used to store motif counts pre- and
    # post-normalization.
    # ({"motif": count(int)}.
    # {"motif": normalized_count(float)}.
    motif_dictionary_pre = {}
    motif_dictionary = {}
    # Used to store nucleotide counts for the file.
    nucleotide_count = 0
    # Opens FASTQ file for reading.
    open_file = open(fastq_file, 'r')
    # Reads in four lines at a time since FASTQ format.
    # Line2 contains nucleotide sequence.
    line1 = open_file.readline().strip("\n")
    line2 = open_file.readline().strip("\n")
    line3 = open_file.readline().strip("\n")
    line4 = open_file.readline().strip("\n")
    # Loops through each line of the FASTQ file. Generates motif
    # counts for each nucleotide line.  
    while len(line2) > 0:
        line_clean = line2.strip("\n")
        # Loops through each nucleotide and generates motif counts. 
        # Generates normalized counts per region and per file.
        for start in range(len(line_clean) - kmer_size + 1):
            kmer = line_clean[start:start + kmer_size]
            current_count = motif_dictionary_pre.get(kmer, 0)
            motif_dictionary_pre[kmer] = current_count + 1
        nucleotide_count = (nucleotide_count + len(line_clean))
        line1 = open_file.readline().strip("\n")
        line2 = open_file.readline().strip("\n")
        line3 = open_file.readline().strip("\n")
        line4 = open_file.readline().strip("\n")
    open_file.close()    
    # Loops through motif_dictionary and normalizes counts 
    # using nucleotide_count from entire file. Also filters
    # out kmers with Ns.   
    for kmer in motif_dictionary_pre:
        if ("N" not in kmer) and ("n" not in kmer):
            normalized_count = (motif_dictionary_pre[kmer] 
                                / float(nucleotide_count))
            motif_dictionary[kmer] = normalized_count
    return motif_dictionary          

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
        file_name = file_split[-1].split(".fastq")[0]
    else:
        file_name = file.split(".fastq")[0]
    return file_name

def normalize_motifs(motif_dictionary, motif_dictionary_sm):
    """Normalizes motif counts using SM input.

        Arguments:
        motif_dictionary -- Dictionary of motif counts for each file.
            {"motif": normalized_count(float)}.
        motif_dictionary_sm -- Dictionary of motif counts for each file.
            {"motif": normalized_count(float)}.

        Output:
        normalized_dictionary -- Dictionary of SM normalized motif counts
            for each file. {"motif": normalized_count(float)}.
    """
    # Used to store normalized motif counts.
    normalized_dictionary = {}
    # Loops through motif_dictionary and normalizes counts using 
    # motif_dictionary_sm. If motif_dictionary_sm does not have
    # motif, skips motif. Skipped motifs should be rare.
    for motif in motif_dictionary:
        if motif_dictionary_sm.get(motif) != None:
            normalized_count = (motif_dictionary[motif] 
                                / float(motif_dictionary_sm[motif]))
            normalized_dictionary[motif] = normalized_count
    return normalized_dictionary


def output_motifs(final_dictionary, file_dictionary, 
                  output_file, file_prefix):
    """Outputs normalized motif enrichment values
        and the file names in header line from
        which they were generated.

        Arguments:
        final_dictionary -- Dictionary of average motif enrichment
            values using all input files. 
            {"motif": normalized_count(float)}.
        file_dictionary -- Dictionary of motif counts for each file.
            {"file_name": {"motif": normalized_count(float)}}.
        output_file-- (-o) Output directory.
        file_prefix -- (-p) Prefix for output files.    

        Output:
        CSV file of RBNS style motif enrichment values for 
            input files. First line contains file names used
            to generate motif enrichment values.
            File name- output_file + file_prefix + "_rbns_motif_enrichment.csv".
            Format - "Motif, Normalized RBNS Motif Enrichment Value".                    
    """
    # Adds "/" to end of output_file if not present.
    if output_file[-1] != "/":
        output_file = output_file + "/" 
    # Used to store output file name.
    output_file_name = (output_file + file_prefix + "_rbns_motif_enrichment.csv")
    # Opens output file for writing.
    output_file = open(output_file_name, 'w')
    # Writes out files used to generate motif enrichment values.
    # Writes header line.
    output_file.write("# Files Used- ")
    for file_name in file_dictionary:
        output_file.write(file_name + "\t")
    output_file.write("\n")    
    output_file.write("Motif, Normalized RBNS Motif Enrichment Value\n")
    # Loops through final_dictionary and outputs motif enrichment values.
    for motif in final_dictionary:
        output_file.write(motif + "," + str(final_dictionary[motif]) + "\n")
    output_file.close()

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
    parser.add_argument("-f", "--fastq", action = "store", type = str, nargs='+', 
                        help="List of paired-end FASTQ files of doseCLIP samples. \
                              Should contain all replicates and read1 and read2 \
                              should be in same order.", 
                        required=True)      
    parser.add_argument("-fs", "--fastq_sm", action = "store", type = str, nargs='+', 
                        help="List of paired-end FASTQ files of SM samples pertaining \
                              to each and all input doseCLIP files. Should be in \
                              same order as doseCLIP samples, including with read1 \
                              and read2.",                        
                        required = True)
    parser.add_argument("-k", "--kmer", action = "store", type = int, 
                        help="Size of kmers to be analyzed.",
                        required = True)
    parser.add_argument("-o", "--output", action = "store", type = str,
                        help="Prefix for output files.",
                        required = True)
    parser.add_argument("-p", "--prefix", action = "store", type = str,
                        help="Prefix for output files.",
                        required = True)                                                                            
    return parser

def main():
    """Main script function. Executes other functions.

        Arguments:
        fastq_list -- (-f) List of paired-end FASTQ files of doseCLIP
            samples. Should contain all replicates and read1 and 
            read2 should be in same order.
        fastq_sm_list -- (-fs) List of paired-end FASTQ files of SM samples 
            pertaining to each and all input doseCLIP files. Should
            be in same order as doseCLIP samples, including with read1
            and read2.
        kmer_size -- (-k) Size of kmer to be analyzed.
        output_file-- (-o) Output directory.
        file_prefix -- (-p) Prefix for output files.

        Output:    
        CSV file of RBNS style motif enrichment values for 
            input files. First line contains file names used
            to generate motif enrichment values.
            File name- output_file + file_prefix + "_rbns_motif_enrichment.csv".
            Format - "Motif, Normalized RBNS Motif Enrichment Value". 
    """
    # Uses argparse to organize command line arguments.
    args_to_parse = init_argparse()
    parsed = args_to_parse.parse_args()
    # Iterates through and executes main functions.
    (final_dictionary, file_dictionary) = get_motif_enrichment(parsed.fastq, 
                                                              parsed.fastq_sm, 
                                                              parsed.kmer)
    output_motifs(final_dictionary, file_dictionary, 
                  parsed.output, parsed.prefix)  

# Executes main function.
if __name__ == "__main__":
    main()