import argparse

def remove_duplicates(blast_out, deseq2_files, output):
    """Removes duplicates using the first row region name
        in the input TSV file and then outputs non-duplicated
        TSV file.

        Arguments:
        blast_out -- List of BLAST produced output files.
        deseq2_files -- List of DESeq2 output files. Used
            to count the number of each type of binding
            region in the output file. Should be in the
            same order as the BLAST input files.
        output -- Output directory.
        
        Output:
        blast_single_region -- Input BLAST files but with only one 
            region per organism comparison.
        counts_file -- Counts of each type of binding region that
            was found in the BLAST output file.
            Format: "Sample,Type,Count\n"     
    """
    # Iterates through each BLAST output file.
    for file in blast_out:
        # Gets BLAST file name.
        if blast_out[0].find("/") != -1:
            file_name = file.split("/")[-1].replace(".tsv", "")
        if output[-1] != "/":
            output += "/"    
        # Opens counts_file and title for output.
        count_file_out = open(output + f"{file_name}_blast_counts_file.tsv", "w")
        count_file_out.write("Sample\tType\tCount\n")
        # sub_gene count dictionary. Used to count the number of 
        # each type of binding region.
        # {"Type":Count}
        sub_gene_count = {}
        # Gets Deseq2 sub-gene data.
        binding_region_dictionary = get_binding_region_to_type(deseq2_files)
        # Set used to check for duplicates.
        blast_single_set = set()
        # Opens BLAST output file.
        with open(file, "r") as blast:
            # Checks for "/" in output directory.
            if output[-1] != "/":
                output += "/"
            # Removes ".tsv" from file name.    
            if file.find(".tsv") != -1:
                file = file.replace(".tsv", "")
            # Creates output file name.
            output_file = (output + file.split("/")[-1].split(".")[0] + "_single_region.tsv")
            # Opens output file.
            with open(output_file, "w") as blast_single_region:
                # Iterates through each line in BLAST output file.
                for line in blast:
                    # Splits line by tab.
                    line = line.strip("\n").split("\t")
                    # If the first line of the file, write to output file.
                    if line[0] not in blast_single_set:
                        blast_single_region.write("\t".join(line) + "\n")
                        blast_single_set.add(line[0])
                        if line[0] in binding_region_dictionary:
                            sub_gene = binding_region_dictionary[line[0]]
                            if sub_gene in sub_gene_count:
                                sub_gene_count[sub_gene] += 1
                            else:
                                sub_gene_count[sub_gene] = 1
                            # Adds count for "all" types.
                            if "all" in sub_gene_count:
                                sub_gene_count["all"] += 1
                            else:
                                sub_gene_count["all"] = 1    
        # Outputs count data.                        
        for sub_gene in sub_gene_count:
            count_file_out.write(f"{output_file}\t{sub_gene}\t{sub_gene_count[sub_gene]}\n")
    # Closes counts file.
    count_file_out.close()        

def get_binding_region_to_type(deseq2_files):
    """Saves each binding region and what genomic annotation
        it is in a dictionary. Used to get counts.

        Arguments:
        deseq2_files -- List of DESeq2 output files. Used
            to count the number of each type of binding
            region in the output file. Should be in the
            same order as the BLAST input files.

        Output:
        binding_region_dictionary -- Dictionary with binding region
            as key and genomic annotation as value.
            {"Region":Type}}        
    """
    # Dictionary with binding region as key and genomic annotation as value.
    binding_region_dictionary = {}
    # Iterates through each DESeq2 output file.
    for file in deseq2_files:
        # Opens DESeq2 output file.
        with open(file, "r") as deseq2:
            # Iterates through each line in DESeq2 output file.
            for line in deseq2:
                # Splits line by tab and saves sub_gene field.
                line_list = line.strip("\n").split(",")
                sub_gene_field = line_list[12]
                # If the first line of the file, skip.
                if line[1] != '"':
                    # Used to save the final sub_gene for output.
                    sub_gene_out = ""
                    # Checks for multiple entries. If so, follows
                    # rules to pick one.
                    if sub_gene_field.find("~") != -1:
                        sub_gene_list = sub_gene_field.split("~")
                        # Search list and finds proper sub-gene output. CDS
                        # overrides all, then exon, and the remaining
                        # would be 5' or 3' UTR.
                        if sub_gene_field.find("3_UTR") != -1:
                            sub_gene_out = "3_UTR"       
                        elif ((sub_gene_field.find("CDS") != -1) or 
                            (sub_gene_field.find("codon") != -1)):
                            sub_gene_out = "exon"
                        elif sub_gene_field.find("5_UTR") != -1:
                            sub_gene_out = "5_UTR"           
                        else:
                            sub_gene_out = "exon"
                    # If no entry, saves as intergenic.     
                    elif len(line_list[7]) > 1:
                        sub_gene_out = "intron"
                    else:
                        sub_gene_out = "intergenic"       
                    # Formats region name to be like BLAST.
                    region_split = line_list[0].replace('"', '').replace("g", "").split("~")
                    region = f"{region_split[0]}:{region_split[1]}-{region_split[2]}({region_split[3]})"
                    # Adds binding region and genomic annotation to dictionary.
                    binding_region_dictionary[region] = sub_gene_out
    return binding_region_dictionary

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
    parser.add_argument("-b", "--blast_out", action = "store", type = str, nargs='+',
                        help="List of BLAST produced output files", required = True)
    parser.add_argument("-d", "--deseq2_files", action = "store", type = str, nargs='+',
                        help="List of DESeq2 output files. Used to count the number of each type of \
                        binding region in the output file. Should be in the same order as the BLAST input files.",
                        required = True)
    parser.add_argument("-o", "--output", action = "store", type = str, default = "",
                        help = "Output directory.", required = True)                                                                          
    return parser

def main():
    """Main script function. Executes other functions.

        Arguments:
        blast_out -- List of BLAST produced output files.
        output -- Output directory. 
        deseq2_files -- List of DESeq2 output files. Used
            to count the number of each type of binding
            region in the output file. Should be in the
            same order as the BLAST input files. 
            
        Output:    
        blast_single_region -- Input BLAST files but with only one 
            region per organism comparison.
        counts_file -- Counts of each type of binding region that
            was found in the BLAST output file.
            Format: "Sample,Type,Count\n"       
    """
    # Uses argparse to organize command line arguments.
    args_to_parse = init_argparse()
    parsed = args_to_parse.parse_args()
    # Iterates through and executes main functions.
    remove_duplicates(parsed.blast_out, parsed.deseq2_files, parsed.output)

# Executes main function.
if __name__ == "__main__":
    main()