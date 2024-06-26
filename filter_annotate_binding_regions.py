import argparse

def filter_annotate_binding_regions(sm_filtered_file, regular_gtf_file, sm_filtered_file_out,
                                    clip_regular_file, clip_regular_file_out):
    """Takes input DeSeq2 produced files (normalized counts or differential
        expression) and executes the data through other functions to output
        the annotated/filtered data.

        Arguments:
        sm_filtered_file -- DeSeq2 produced SM counts file (CSV) for a single protein
            concentration (produced with replicates). This should have been analyzed
            and normalized with DeSeq2 with the regular CLIP samples and the the SM 
            Input samples. Also, this file could be the filtered differential expression
            output file or the normalized counts file.
        regular_gtf_file -- Gencode/Ensembl/UCSC provided GTF file. Should be the GTF
            annotation for the target alignment sequence. The file is used to annotate
            the binding regions.
        sm_filtered_file_out -- Output normalized SM filtered counts or SM differential 
            expression file (this is output from sm_filtered_file). The file will 
            be annotated. The file is in CSV format.    
        clip_regular_file -- DeSeq2 produced CLIP counts file (CSV) for a single protein
            concentration (produced with replicates). This file should have been analyzed
            and normalized with DeSeq2 without the SM Input samples. Also, this file could 
            be the filtered differential expression output file or the normalized
            counts file. These are the numerical values that are primarily used for later
            analyses. This is an OPTIONAL parameter, as the script is able to just annotate 
            the SM filtered counts as well.
        clip_regular_file_out -- Output normalized regular CLIP filtered counts or 
            CLIP differential expression file (this is output from clip_regular_file). 
            The file will be annotated and filtered with the SM filtered binding regions.
            The file is in CSV format. This is an OPTIONAL parameter, as the script is 
            able to just annotate the SM filtered counts as well.

        Output:
            This function has no output.        

    """
    # Processes SM filtered file for output.
    output_file(sm_filtered_file, regular_gtf_file, sm_filtered_file_out,
                "SM", clip_regular_file, clip_regular_file_out)
    # Checks to see if clip_regular_file was input. If so, processes it
    # for filtering and annotation.   
    if clip_regular_file != None:
        output_file(sm_filtered_file, regular_gtf_file, sm_filtered_file_out,
                     "CLIP", clip_regular_file, clip_regular_file_out)   
                        
def add_events_to_dictionary(filtered_file, regular_gtf_file):
    """Takes an input counts file, annotates it and adds to a dictionary
        for output.

        Arguments:
        filtered_file -- DeSeq2 produced SM counts file or CLIP file (CSV) for a 
            single protein concentration (produced with replicates). This should have 
            been analyzed and normalized with DeSeq2 with the regular CLIP samples and 
            the the SM Input samples. Also, this file could be the filtered 
            differential expression output file or the normalized counts file.
        regular_gtf_file -- Gencode/Ensembl/UCSC provided GTF file. Should be the GTF
            annotation for the target alignment sequence. The file is used to annotate
            the binding regions.
            
        Output:    
        region_dictionary -- Filled output dictionary with each binding region names
            as the key and the clean line as the value.
            {"region_name": "line_clean=gene_line=sub_gene_line"}.
        title_line -- String containing the title line of file. Used for output.

    """                      
    # Dictionary that stores each read pileup region.
    region_dictionary = {}
    # Saves title line.
    title_line = ""
    # Saves GTF information to object for faster annotation.
    gtf_object = make_gtf_object(regular_gtf_file)
    # Iterates through each line in the counts file (CSV) and adds it to the dictionary.
    # Also, looks for gene and sub-gene information from GTF file.                     
    with open(filtered_file) as open_filtered_file:
        for line in open_filtered_file:
            # Cleans and parses line.
            line_clean = line.strip("\n")
            line_split = line_clean.split(",")
            # Looks for and saves title line. Passes line after saving.
            if (len(line_split[0]) == 0 or (line_clean[1] == '"')):
                title_line = line_clean
                continue
            # Used to save gene and sub-gene information, for addition
            # to output after iterating through GTF file.
            gene_line = ""
            sub_gene_line = ""
            gene_list = []
            # Parses gene information.
            gene_split = line_split[0].split("~")
            # Saves information needed to compare to GTF.
            chromosome, cordinate1 = gene_split[0].replace('"',''), int(gene_split[1])
            cordinate2, strand = int(gene_split[2]), gene_split[3][0]
            chrom_dictionary = gtf_object.get(chromosome)
            if chrom_dictionary != None:
                gene_list = chrom_dictionary.get(strand)
            if gene_list != None:
                for gene in gene_list:
                    if (((gene[0] <= cordinate1) and (cordinate1 <= gene[1])) or
                       ((gene[0] <= cordinate2) and (cordinate2 <= gene[1])) or
                       ((gene[0] <= cordinate1) and (cordinate2 <= gene[1])) or
                       ((gene[0] >= cordinate1) and (cordinate2 >= gene[1]))):
                        gene_id, gene_type, gene_name, event_type = gene[2], gene[3], gene[4], gene[5]
                        # Checks to see if there has been a gene that has already been
                        # found that overlaps. This should be a rare occurance. If no
                        # previous gene, saves new data, if previous then adds additional
                        # gene information to string for output.
                        if (event_type == "gene"):
                            if len(gene_line) == 0:
                                gene_line = f"{gene_id}~{gene_type}~{gene_name}"
                            else:
                                gene_line = f"{gene_line}~{gene_id}~{gene_type}~{gene_name}"   
                        # Performs same as above for the gene line but for sub-gene features.            
                        if ((event_type != "transcript") and (event_type != "gene")):
                        # Classifies UTR type. If first coordinate is closer to first transcript
                        # coordinate, classifies as a 5' UTR, else classifies as 3' UTR.     
                            if len(sub_gene_line) == 0:
                                sub_gene_line = event_type
                            else:    
                                sub_gene_line = f"{sub_gene_line}~{event_type}"
            # Adds all information to region_dictionary for output.
            region_dictionary[line_split[0]] = f"{line_clean}={gene_line}={sub_gene_line}"         
    open_filtered_file.close()
    return region_dictionary, title_line

def output_file(sm_filtered_file, regular_gtf_file, sm_filtered_file_out, sample_keyword,
                clip_regular_file = "None", clip_regular_file_out = "None"):
    """Takes DeSeq2 produced files, annotates, filters if SM and regular
        CLIP files are input, and outputs in CSV format. If only SM 
        filtered file is input, only annotates and outputs. 

        Arguments:
        sm_filtered_file -- DeSeq2 produced SM counts file (CSV) for a single protein
            concentration (produced with replicates). This should have been analyzed
            and normalized with DeSeq2 with the regular CLIP samples and the the SM 
            Input samples. Also, this file could be the filtered differential expression
            output file or the normalized counts file. 
        clip_regular_file -- DeSeq2 produced CLIP counts file (CSV) for a single protein
            concentration (produced with replicates). This file should have been analyzed
            and normalized with DeSeq2 without the SM Input samples. Also, this file could 
            be the filtered differential expression output file or the normalized
            counts file. These are the numerical values that are primarily used for later
            analyses. This is an OPTIONAL parameter, as the script is able to just annotate 
            the SM filtered counts as well.
        regular_gtf_file -- Gencode/Ensembl/UCSC provided GTF file. Should be the GTF
            annotation for the target alignment sequence. The file is used to annotate
            the binding regions.
        sample_keyword -- Keyword to tell the function to process the sm_filtered_file or
            the clip_regular_file. Keyword options are "SM" for sm_filtered_file and
            "CLIP" for clip_regular_file.    
            
        Output:    
        sm_filtered_file_out -- Output normalized SM filtered counts or SM differential 
            expression file (this is output from sm_filtered_file). The file will 
            be annotated. The file is in CSV format. 
            File format (Title line is orignal from file):
            title_line, gene_id, gene_type, gene_name, sub_gene_type, all_genes, all_sub_gene_types
        clip_regular_file_out -- Output normalized regular CLIP filtered counts or 
            CLIP differential expression file (this is output from clip_regular_file). 
            The file will be annotated and filtered with the SM filtered binding regions.
            The file is in CSV format. This is an OPTIONAL parameter, as the script is 
            able to just annotate the SM filtered counts as well.
            File format (Title line is orignal from file):
            title_line, gene_id, gene_type, gene_name, sub_gene_type, all_genes, all_sub_gene_types
    """
    # Initializes variables needed to be iterated through.
    region_dictionary_sm = {}
    region_dictionary = {}
    # If CLIP data dictionary is ued for filtered data, if SM outputs all events.
    region_dictionary_count = {}
    # Stores annotated binding regions in dictionary and writes out the modified title line.
    if sample_keyword == "SM":
        region_dictionary_count, title_line = add_events_to_dictionary(sm_filtered_file, regular_gtf_file)
        # Modifies title line for new annotated output.
        title_line = f"{title_line},gene_id,gene_type,gene_name,sub_gene_type,all_genes,all_sub_gene_types"
        file_out_open = open(sm_filtered_file_out, "w")
        file_out_open.write(title_line + "\n")
        # Counts and outputs sub-gene types for the sample.
        count_and_output_sub_genes(region_dictionary_count, sm_filtered_file_out, None)
    # Stores annotated binding regions in dictionary and writes out the modified title line.    
    elif sample_keyword == "CLIP":
        region_dictionary, title_line = add_events_to_dictionary(clip_regular_file, regular_gtf_file)
        # Also grabs SM filtered events so the dictionary can be used for filtering output.
        region_dictionary_sm, title_line = add_events_to_dictionary(sm_filtered_file, regular_gtf_file)
        # Modifies title line for new annotated output.
        title_line = f"{title_line},gene_id,gene_type,gene_name,sub_gene_type,all_genes,all_sub_gene_types"
        file_out_open = open(clip_regular_file_out, "w")
        file_out_open.write(title_line + "\n")
        # Filters CLIP dictionary with SM data.
        for event in region_dictionary:
            if event in region_dictionary_sm:
                region_dictionary_count[event] = (region_dictionary.get(event))
        # Counts and outputs sub-gene types for the sample.
        count_and_output_sub_genes(region_dictionary_count, None, clip_regular_file_out)
    # Iterates through dictionary and outputs each annotated binding region.           
    for region in region_dictionary_count:
        # Strings used to output multiple genes or sub-genes. 
        # This should be rare.
        all_genes = ""
        all_sub_genes = ""
        region_line = region_dictionary_count.get(region)
        region_split = region_line.split("=")
        # Checks for multiple gene entries and outputs appropriately.
        if (region_split[1].find("~") != -1):
            gene_split = region_split[1].split("~")
            gene = f"{gene_split[0]},{gene_split[1]},{gene_split[2]}"
            if len(gene_split) > 3:
                all_genes = region_split[1]
        else:
            # Checks for no annotated genes. Adds commas so columns will
            # be an empty value if if no genes were found.
            if len(region_split[1]) != 0:
                gene = region_split[1]
            else:
                gene = ",,"
        # Checks for multiple sub-gene entries and outputs appropriately.
        if (region_split[2].find("~") != -1):
            sub_gene_split = region_split[2].split("~")
            sub_gene = f"{sub_gene_split[0]}"
            all_sub_genes = region_split[2]
        else:               
            sub_gene = region_split[2]
        output = f"{region_split[0]},{gene},{sub_gene},{all_genes},{all_sub_genes}"
        file_out_open.write(output + "\n")
    file_out_open.close()

def count_and_output_sub_genes(region_dictionary, sm_filtered_file_out, clip_regular_file_out):
    """Takes an input region dictionary and outputs counts of annotation information.

        Arguments:
        region_dictionary -- Filled output dictionary with each binding region names
            as the key and the clean line as the value.
            {"region_name": "line_clean=gene_line=sub_gene_line"}.
            
        Output:    
        sm_filtered_file_out -- Output of annotation counts for normalized SM filtered 
            counts or SM differential expression file (this is output from sm_filtered_file). 
            The file is in CSV format. File contains a ".annotation_counts.csv" suffix. If no
            input (both files will not be input at once), None is input and output will be
            ignored.
            File Format: Sub-Gene Type, Sub-Gene Count
        clip_regular_file_out -- Output of annotation counts for normalized regular CLIP 
            filtered counts or CLIP differential expression file (this is output from 
            clip_regular_file). The file will be annotated and filtered with the SM 
            filtered binding regions. The file is in CSV format. This is an OPTIONAL 
            parameter, as the script is able to just annotate the SM filtered counts 
            as well. File contains a ".annotation_counts.csv" suffix. If no input 
            (both files will not be input at once), None is input and output will be
            ignored.
            File Format: Sub-Gene Type, Sub-Gene Count
    """
    # Initializes dictionary for sub-gene regions and counts.
    # [sub_gene] = sub_gene_count
    sub_gene_dictionary = {}
    # Checks for the file that was input. Adds the ".annotation_counts.csv" suffix to the
    # output file name and opens the file for writing.
    if clip_regular_file_out == None:
        sm_filtered_file_out = sm_filtered_file_out.replace(".csv", ".annotation_counts.csv")
        file_out_open = open(sm_filtered_file_out, "w")
    else:
        clip_regular_file_out = clip_regular_file_out.replace(".csv", ".annotation_counts.csv")
        file_out_open = open(clip_regular_file_out, "w")
    # Iterates through dictionary and outputs each annotated binding region.          
    for region in region_dictionary:
        # Sets gene_flag boolean to indicate if region was present in
        # a gene.
        gene_flag = False
        # Strings used to output multiple genes or sub-genes. 
        # This should be rare.
        all_genes = ""
        all_sub_genes = ""
        region_line = region_dictionary.get(region)
        region_split = region_line.split("=")
        # Checks for multiple gene entries and saves as variables. Not all
        # variables are output at this point. Saved for later script 
        # modifications.
        if (region_split[1].find("~") != -1):
            gene_flag = True
            gene_split = region_split[1].split("~")
            gene = f"{gene_split[0]},{gene_split[1]},{gene_split[2]}"    
            if len(gene_split) > 3:
                all_genes = region_split[1]
        else:
            # Checks for no annotated genes. Adds commas so columns will
            # be an empty value if if no genes were found.
            # Sets boolean gene_flag to false because no gene was annotated
            # for the binding region.
            if len(region_split[1]) != 0:
                gene_flag = True
                gene = region_split[1]
            else:
                gene = ",,"
        # Checks sub-genes and categorizes and counts for proper output.
        # Checks for multiple sub-gene entries and iterates through, using
        # rules to classify and catagorize the sub-gene appropriately.
        if (region_split[2].find("~") != -1):
            sub_gene_split = region_split[2].split("~")
            for sub_gene_list in sub_gene_split:
                # If CDS, then overrides other sub-types.
                if sub_gene_list == "CDS":
                    sub_gene = "CDS"
                    break
                # If start or stop codon, this is just classified as an exon.
                elif ((sub_gene_list == "start_codon") or (sub_gene_list == "stop_codon")):
                    sub_gene = "exon"
                # Only remaining should be exon.
                else:
                   sub_gene = sub_gene_list      
            # Saves all in-case of future code revisions and use.       
            all_sub_genes = region_split[2]
        else:               
            sub_gene = region_split[2]
        # If in a gene but no sub-gene type, classified as an intron.    
        if ((len(sub_gene) == 0) and (gene_flag == True)):
            sub_gene = "intron"
        # If not in gene or intron, classified as intergenic.     
        elif ((len(sub_gene) == 0) and (gene_flag == False)):
            sub_gene = "intergenic"
        # Adds sub-gene to dictionary with proper count as key.
        if sub_gene not in sub_gene_dictionary:
             sub_gene_dictionary[sub_gene] = 1
        else:
             sub_gene_count = sub_gene_dictionary.get(sub_gene)
             sub_gene_count = (sub_gene_count + 1)
             sub_gene_dictionary[sub_gene] = sub_gene_count       
    # Writes out each sub-gene type and count.         
    file_out_open.write("Sub-Gene Type,Sub-Gene Count\n")
    for sub_gene in sub_gene_dictionary:
        sub_gene_count = sub_gene_dictionary.get(sub_gene)
        output = f"{sub_gene},{sub_gene_count}\n"
        file_out_open.write(output)
    file_out_open.close()

def make_gtf_object(regular_gtf_file):
    """Takes a GTF file and creates a GTF object from it.

        Arguments:
        regular_gtf_file -- Gencode/Ensembl/UCSC provided GTF file. Should be the GTF
            annotation for the target alignment sequence. The file is used to annotate
            the binding regions.
            
        Output:    
        gtf_object -- GTF object created from the GTF file. Object is organized by
            chromosome, strand, and then cordinate. The object is used to annotate
            the binding regions.
        {chromsome: {strand: [coordinate1, coordinate2, gene_id, gene_type, 
                              gene_name, event_type]}}
    """
    gtf_object = {}
    # Iterates through each line in the GTF file (TSV) and adds 
    # each line to the GTF object.                      
    with open(regular_gtf_file) as open_regular_gtf_file:
        for line_gtf in open_regular_gtf_file:
            # Skips title lines.
            if line_gtf[0] != "#":
                line_gtf_split = line_gtf.strip("\n").split("\t")
                # Saves data needed for comparison.
                chromosome, coordinate1 = line_gtf_split[0], int(line_gtf_split[3])
                coordinate2, strand = int(line_gtf_split[4]), line_gtf_split[6]
                # Saves the event feature (type) as a variable.
                event_type = line_gtf_split[2]
                # If line is a gene line, saves needed gene information to string for output. 
                if (event_type == "gene"):
                    gene_data = line_gtf_split[8].split(";")
                    gene_id_pre = gene_data[0].split('"')
                    gene_id = gene_id_pre[1]
                    gene_type_pre = gene_data[1].split('"')
                    gene_type = gene_type_pre[1]
                    gene_name_pre = gene_data[2].split('"')
                    gene_name = gene_name_pre[1]
                # Performs same as above for the gene line but for sub-gene features.            
                if ((event_type != "transcript") and (event_type != "gene")):
                    # Classifies UTR type. If first coordinate is closer to first transcript
                    # coordinate, classifies as a 5' UTR, else classifies as 3' UTR.
                    if (event_type == "UTR"):
                        utr_cord1 = int(line_gtf_split[3])
                        utr_cord2 = int(line_gtf_split[4])
                        if strand == "+":
                            cord1_difference = abs(transcript_cord1 - utr_cord1)
                            cord2_difference = abs(transcript_cord2 - utr_cord1)
                            if (cord1_difference > cord2_difference):
                                event_type = "3_UTR"
                            else:
                                event_type = "5_UTR"
                        else:
                            cord1_difference = abs(transcript_cord1 - utr_cord2)
                            cord2_difference = abs(transcript_cord2 - utr_cord2)
                            if (cord1_difference > cord2_difference):
                                event_type = "5_UTR"
                            else:
                                event_type = "3_UTR"        
                    # If transcript feature, saves coordinates to determine UTR type.        
                if (event_type == "transcript"):
                    transcript_cord1 = int(line_gtf_split[3])
                    transcript_cord2 = int(line_gtf_split[4])             
                # Adds all information to gtf_object for output. Has to check to see if
                # keys have been added previously. If so, adds to existing key, if not
                # creates new key.
                if chromosome in gtf_object:
                    if strand in gtf_object[chromosome]:
                        gtf_object[chromosome][strand].append([coordinate1, coordinate2, gene_id, 
                                                               gene_type, gene_name, event_type])
                    else:
                        gtf_object[chromosome][strand] = [[coordinate1, coordinate2, gene_id, 
                                                           gene_type, gene_name, event_type]]
                else:
                    gtf_object[chromosome] = {strand: [[coordinate1, coordinate2, gene_id,
                                                        gene_type, gene_name, event_type]]}
    open_regular_gtf_file.close()
    return gtf_object            

def init_argparse():
    """Initiates the use of argparse. Returns parsed
        command line arguments.

        Arguments:
        See command line argument descriptions below.

        Output:
        Returns parser object. Object is used for function input.
    """
    # Script description.
    parser = argparse.ArgumentParser(description = "Takes input DeSeq2 produced files (normalized counts \
                                                    or differential expression) and annotates the files, \
                                                    and if normal CLIP file was input, filters it using \
                                                    the SM Input file.")        
    # Command line arguments and descriptions.
    parser.add_argument("-sm", "--sm_filtered_file", action = "store", type = str, 
                        help="DeSeq2 produced SM counts file (CSV) for a single protein \
                              concentration (produced with replicates). This should \
                              have been analyzed and normalized with DeSeq2 with the \
                              regular CLIP samples and the the SM Input samples. Also, \
                              this file could be the filtered differential expression \
                              output file or the normalized counts file.", required = True)
    parser.add_argument("-gtf", "--regular_gtf_file", action = "store", type = str, 
                        help = "Gencode/Ensembl/UCSC provided GTF file. Should be the GTF \
                                annotation for the target alignment sequence. The file is \
                                used to annotate the binding regions.", required = True)
    parser.add_argument("-clip", "--clip_regular_file", action = "store", type = str,
                        help = "DeSeq2 produced CLIP counts file (CSV) for a single protein \
                                concentration (produced with replicates). This file should \
                                have been analyzed and normalized with DeSeq2 without the \
                                SM Input samples. Also, this file could be the filtered \
                                differential expression output file or the normalized counts \
                                file. These are the numerical values that are primarily used \
                                for later analyses. This is an OPTIONAL parameter, as the \
                                script is able to just annotate the SM filtered counts as well.", 
                                required = False)
    parser.add_argument("-o_sm", "--sm_filtered_file_out", action = "store", type = str, 
                        help="Output annotated SM file in CSV format. This is output \
                              from the sm_filtered_file input.", 
                        required = True)
    parser.add_argument("-o_cl", "--clip_regular_file_out", action = "store", type = str, 
                        help="Output annotated regular CLIP file in CSV format. This is \
                              output from the clip_regular_file input.", 
                        required = False)                                                                              
    return parser

def main():
    """Main script function. Executes other functions.

        Arguments:
        sm_filtered_file -- DeSeq2 produced SM counts file (CSV) for a single protein
            concentration (produced with replicates). This should have been analyzed
            and normalized with DeSeq2 with the regular CLIP samples and the the SM 
            Input samples. Also, this file could be the filtered differential expression
            output file or the normalized counts file.
        regular_gtf_file -- Gencode/Ensembl/UCSC provided GTF file. Should be the GTF
            annotation for the target alignment sequence. The file is used to annotate
            the binding regions.    
        clip_regular_file -- DeSeq2 produced CLIP counts file (CSV) for a single protein
            concentration (produced with replicates). This file should have been analyzed
            and normalized with DeSeq2 without the SM Input samples. Also, this file could 
            be the filtered differential expression output file or the normalized
            counts file. These are the numerical values that are primarily used for later
            analyses. This is an OPTIONAL parameter, as the script is able to just annotate 
            the SM filtered counts as well.
            
        Output:    
        sm_filtered_file_out -- Output normalized SM filtered counts or SM differential 
            expression file (this is output from sm_filtered_file). The file will be annotated.
            The file is in CSV format.
        clip_regular_file_out -- Output normalized regular CLIP filtered counts or CLIP differential 
            expression file (this is output from clip_regular_file). The file will be annotated and
            filtered with the SM filtered binding regions. The file is in CSV format.
            This is an OPTIONAL parameter, as the script is able to just annotate 
            the SM filtered counts as well.    
    """
    # Uses argparse to organize command line arguments.
    args_to_parse = init_argparse()
    parsed = args_to_parse.parse_args()
    # Iterates through and executes main functions.
    filter_annotate_binding_regions(parsed.sm_filtered_file, parsed.regular_gtf_file,
                                    parsed.sm_filtered_file_out, parsed.clip_regular_file,
                                    parsed.clip_regular_file_out)

# Executes main function.
if __name__ == "__main__":
    main()