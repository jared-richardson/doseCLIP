import argparse        
                        
def add_events_to_dictionary(bed_file, regular_gtf_file):
    """Takes input BED file with binding regions and annotated each region using
        a GTF file. Outputs information in dictionary.

        Arguments:
        bed_file -- BED file with binding regions. This file
            would be produced from a separate experiment where there
            where no SM Input samples.
        regular_gtf_file -- Gencode/Ensembl/UCSC provided GTF file. Should be the GTF
            annotation for the target alignment sequence. The file is used to annotate
            the binding regions.

        Output:
        region_dictionary -- Dictionary containing the event name (chromosome information)
            as the key and the entire BED line with gene information as the value.
            {"chromosome_gtf,cordinate1_gtf,cordinate2_gtf,strand_gtf": 
             "line_clean=gene_line=sub_gene_line"}   
    """
    # Dictionary that stores each read pileup region.
    region_dictionary = {}
    # Iterates through the BED file and saves chromosome information for comparison
    # to GTF file. Compares each event to each event in the GTF file and adds gene
    # information with the event information to a dictionary for output.                  
    with open(bed_file) as open_bed_file:
        for line in open_bed_file:
            # Cleans and parses line.
            line_clean = line.strip("\n").replace("\t",",")
            line_split = line_clean.split(",")
            # Used to save gene and sub-gene information, for addition
            # to output after iterating through GTF file.
            gene_line = ""
            sub_gene_line = ""
            # Saves information needed to compare to GTF.
            chromosome, cordinate1, cordinate2, strand = (line_split[0], int(line_split[1]),
                                                          int(line_split[2]), line_split[5])
            event_name = f"{chromosome}~{cordinate1}~{cordinate2}~{strand}" 
            # Iterates through each line in the GTF file (TSV) and looks for overlap
            # with binding region. If overlap, saves needed information for output.                      
            with open(regular_gtf_file) as open_regular_gtf_file:
                for line_gtf in open_regular_gtf_file:
                    if line_gtf[0] != "#":
                        line_gtf_split = line_gtf.strip("\n").split("\t")
                        # Saves data needed for comparison.
                        chromosome_gtf, cordinate1_gtf = line_gtf_split[0], int(line_gtf_split[3])
                        cordinate2_gtf, strand_gtf = int(line_gtf_split[4]), line_gtf_split[6]
                        # Checks for identical chromosome, strand, and overlap between the
                        # binding region and the gene/sub-gene feature. 
                        if ((chromosome == chromosome_gtf) and (strand == strand_gtf)
                            and (((cordinate1_gtf <= cordinate1) and (cordinate1 <= cordinate2_gtf)) 
                                or ((cordinate1_gtf <= cordinate2) and (cordinate2 <= cordinate2_gtf))
                                or ((cordinate1_gtf <= cordinate1) and (cordinate2 <= cordinate2_gtf))
                                or ((cordinate1_gtf >= cordinate1) and (cordinate2 >= cordinate2_gtf)))):
                            # If GTF is a gene line, saves needed gene information to string for output. 
                            if (line_gtf_split[2] == "gene"):
                                gene_data = line_gtf_split[8].split(";")
                                gene_id_pre = gene_data[0].split('"')
                                gene_id = gene_id_pre[1]
                                gene_type_pre = gene_data[1].split('"')
                                gene_type = gene_type_pre[1]
                                gene_name_pre = gene_data[2].split('"')
                                gene_name = gene_name_pre[1]
                                # Checks to see if there has been a gene that has already been
                                # found that overlaps. This should be a rare occurance. If no
                                # previous gene, saves new data, if previous then adds additional
                                # gene information to string for output.
                                if len(gene_line) == 0:
                                    gene_line = f"{gene_id}~{gene_type}~{gene_name}"
                                else:
                                    gene_line = f"{gene_line}~{gene_id}~{gene_type}~{gene_name}"
                            # Performs same as above for the gene line but for sub-gene features.            
                            elif (line_gtf_split[2] != "transcript"):
                                if len(sub_gene_line) == 0:
                                    sub_gene_line = line_gtf_split[2]
                                else:    
                                    sub_gene_line = f"{sub_gene_line}~{line_gtf_split[2]}"
            # Adds all information to region_dictionary for output.
            region_dictionary[event_name] = f"{line_clean}={gene_line}={sub_gene_line}" 
    open_bed_file.close()
    return region_dictionary

def annotate_output_file(bed_file, regular_gtf_file, bed_file_out):
    """Takes input BED file with binding regions and annotated each region using
        a GTF file. Outputs information in BED format with additional columns.

        Arguments:
        bed_file -- BED file with binding regions. This file
            would be produced from a separate experiment where there
            where no SM Input samples.
        regular_gtf_file -- Gencode/Ensembl/UCSC provided GTF file. Should be the GTF
            annotation for the target alignment sequence. The file is used to annotate
            the binding regions.

        Output:
        bed_file_out -- GTF Annotated BED file with binding regions. This file
            would be produced from a separate experiment where there
            where no SM Input samples.
    """
    # Opens outfile for writing to.
    bed_file_out_write = open(bed_file_out, "w")
    # Finds overlapping gene events and adds them to a dictionary.
    region_dictionary = add_events_to_dictionary(bed_file, regular_gtf_file) 
    # Counts and outputs sub-gene types for the sample.
    count_and_output_sub_genes(region_dictionary, bed_file_out) 
    # Iterates through dictionary and outputs each annotated binding region.           
    for region in region_dictionary:
        # Strings used to output multiple genes or sub-genes. 
        # This should be rare.
        all_genes = ""
        all_sub_genes = ""
        region_line = region_dictionary.get(region)
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
        # Adds tabs for output.
        output_tab = output.replace(",","\t")
        bed_file_out_write.write(output_tab + "\n")
    bed_file_out_write.close()

def count_and_output_sub_genes(region_dictionary, bed_file_out):
    """Takes an input region dictionary and outputs counts of annotation information.

        Arguments:
        region_dictionary -- Filled output dictionary with each binding region names
            as the key and the clean line as the value.
            {"region_name": "line_clean=gene_line=sub_gene_line"}.
            
        Output:    
        bed_file_out -- Sub-gene counts of binding regions. This file
            would be produced from a separate experiment where there
            where no SM Input samples.
            File Format: Sub-Gene Type, Sub-Gene Count
    """
    # Initializes dictionary for sub-gene regions and counts.
    # [sub_gene] = sub_gene_count
    sub_gene_dictionary = {}
    # Adds the ".annotation_counts.csv" suffix to the
    # output file name and opens the file for writing.
    bed_file_out = bed_file_out.replace(".bed", ".annotation_counts.csv")
    file_out_open = open(bed_file_out, "w")
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

def init_argparse():
    """Initiates the use of argparse. Returns parsed
        command line arguments.

        Arguments:
        See command line argument descriptions below.

        Output:
        Returns parser object. Object is used for function input.
    """
    # Script description.
    parser = argparse.ArgumentParser(description = "Takes input BED file with binding regions and \
                                                    annotated each region using a GTF file. Outputs \
                                                    information in BED format with additional columns.")        
    # Command line arguments and descriptions.
    parser.add_argument("-bed", "--bed_file", action = "store", type = str, 
                        help="BED file with binding regions. This file \
                             would be produced from a separate experiment where there \
                             where no SM Input samples..", required = True)
    parser.add_argument("-gtf", "--regular_gtf_file", action = "store", type = str, 
                        help = "Gencode/Ensembl/UCSC provided GTF file. Should be the GTF \
                                annotation for the target alignment sequence. The file is \
                                used to annotate the binding regions.", required = True)
    parser.add_argument("-o", "--bed_file_out", action = "store", type = str, 
                        help="GTF Annotated BED file with binding regions. This file \
                              would be produced from a separate experiment where there \
                              where no SM Input samples.", 
                        required = True)                                                                            
    return parser

def main():
    """Main script function. Executes other functions.

        Arguments:
        bed_file -- BED file with binding regions. This file
            would be produced from a separate experiment where there
            where no SM Input samples.
        regular_gtf_file -- Gencode/Ensembl/UCSC provided GTF file. Should be the GTF
            annotation for the target alignment sequence. The file is used to annotate
            the binding regions.

        Output:
        bed_file_out -- GTF Annotated BED file with binding regions. This file
            would be produced from a separate experiment where there
            where no SM Input samples.    
    """
    # Uses argparse to organize command line arguments.
    args_to_parse = init_argparse()
    parsed = args_to_parse.parse_args()
    # Iterates through and executes main functions.
    annotate_output_file(parsed.bed_file, parsed.regular_gtf_file,
                         parsed.bed_file_out)

# Executes main function.
if __name__ == "__main__":
    main()