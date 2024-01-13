"""test_motif_secondary_features.py- Unit tests for
    motif_secondary_features.py
"""
import pytest
# Imports all functions from make_random_bed_for_motif.py.
import motif_secondary_features_structure as m_s_f_s

"""
Tests organize_nucleotides with various nucleotide lists.

    Data Types:
    fasta_file_name -- String, name of FASTA file being analyzed
    kmer_size -- Size of kmer to be analyzed.
    motif_list -- List of motifs to be analyzed. Default is YGCY.
    motif_title -- Title for motif output file. Default is YGCY.
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
                         region_size]
    structure_dictionary -- Dictionary of secondary structure predictions for
        motifs and surrounding nucleotides. The data is a dictionary with the
        region as the key and the list of secondary structure predictions as
        the value, arranged in three lists.
        {"region_name": [structure, y_structure, gc_structure]}. 
"""

@pytest.mark.parametrize("fasta_file_name, kmer_size, motif_list, \
                          motif_title, region_dictionary, structure_dictionary", [
                        # Test 1: Single region with two motif groups. One
                        # group with two motifs, the other with one. Also
                        # tests a motif coming up first and last in the
                        # sequence.
                        # fasta_file_name
                        (("test_files/motif_secondary_features_structure/test1.fasta"),
                        # kmer_size
                        (4),
                        # motif_list
                        (["CGCC", "CGCU", "UGCC", "UGCU"]),
                        # motif_title
                        ("YGCY"),
                        # region_dictionary
                        ({"name1": [['AAAAUUUUUUUUUUUUUU'], 
                                    ['AAAA'], ['UUUUUUUU'], 
                                    ['UGCUAAAACGCCUUUUUUU', 'UUUUUUUCGCC'], 
                                    24]}),
                        # structure_dictionary
                        ({"name1": [['(','(','(','(','(','(','(','(','(','(','(','(','(','(','(','(','(','.'], 
                                    ['(','(','(','(','.','.'], ['(','(','(','(','.','.']]})),
                        # Test 2: Single region with two motif groups. One
                        # group with two motifs, the other with a double
                        # motif. Also checks for motifs that are not at
                        # beginning or end but the surrounding nucleotides
                        # would be.
                        # fasta_file_name
                        (("test_files/motif_secondary_features_structure/test2.fasta"),
                        # kmer_size
                        (4),
                        # motif_list
                        (["CGCC", "CGCU", "UGCC", "UGCU"]),
                        # motif_title
                        ("YGCY"),
                        # region_dictionary
                        ({"name1": [['CCCAAAAUUUUUUUUUUUUUUGGGG'], 
                                    ['AAAA'], ['UUUUUUUU'], 
                                    ['CCCUGCUAAAACGCCUUUUUUU', 
                                     'UUUUUUUCGCCCGCCGGGG'], 
                                    35]}),
                        # structure_dictionary
                        ({"name1": [['.','.','.',')',')',')',')','(','(','(','(','(','(','(','(','(','(','(','(','(','(',')',')',')',')'], 
                                    ['(',')','.','.',')',')','(',')'], ['(',')','.','.',')',')','.','.']]})),
                        # Test 3: Multiple region test. This is just 
                        # the regions from Test 1 and Test 2 combined.
                        # fasta_file_name
                        (("test_files/motif_secondary_features_structure/test3.fasta"),
                        # kmer_size
                        (4),
                        # motif_list
                        (["CGCC", "CGCU", "UGCC", "UGCU"]),
                        # motif_title
                        ("YGCY"),
                        # region_dictionary
                        ({"name1": [['AAAAUUUUUUUUUUUUUU'], 
                                    ['AAAA'], ['UUUUUUUU'], 
                                    ['UGCUAAAACGCCUUUUUUU', 'UUUUUUUCGCC'], 
                                    24],
                          "name2": [['CCCAAAAUUUUUUUUUUUUUUGGGG'], 
                                    ['AAAA'], ['UUUUUUUU'], 
                                    ['CCCUGCUAAAACGCCUUUUUUU', 
                                     'UUUUUUUCGCCCGCCGGGG'], 
                                    35]}),
                        # structure_dictionary
                        ({"name1": [['(','(','(','(','(','(','(','(','(','(','(','(','(','(','(','(','(','.'], 
                                    ['(','(','(','(','.','.'], ['(','(','(','(','.','.']],
                          "name2": [['.','.','.',')',')',')',')','(','(','(','(','(','(','(','(','(','(','(','(','(','(',')',')',')',')'], 
                                    ['(',')','.','.',')',')','(',')'], ['(',')','.','.',')',')','.','.']]})),
                        # Test 4: Multiple lines with no motifs.
                        # fasta_file_name
                        (("test_files/motif_secondary_features_structure/test4.fasta"),
                        # kmer_size
                        (4),
                        # motif_list
                        (["CGCC", "CGCU", "UGCC", "UGCU"]),
                        # motif_title
                        ("YGCY"),
                        # region_dictionary
                        ({"name1": [[], [], [], [], 22],
                          "name6": [[], [], [], [], 9]}),
                        # structure_dictionary
                        ({"name1": [[], [], []],
                          "name6": [[], [], []]})),                                                                      
                        ])

def test_organize_nucleotides(fasta_file_name, kmer_size, motif_list,
                              motif_title, region_dictionary, 
                              structure_dictionary):
    """
        GIVEN a list of FASTA files with secondary structure predictions.
        WHEN the function extracts the secondary structure information
            for the motifs and around the motifs
        THEN the output dictionary is checked for the correct data.
    """
    # Executes function and checks output.
    assert m_s_f_s.organize_nucleotides(fasta_file_name, kmer_size, motif_list,
                                      motif_title) == (region_dictionary, 
                                                       motif_title, motif_list,
                                                       structure_dictionary)


"""
Tests output_region_and_file_data with various seondary structure dictionaries.

    Data Types:
    file_dictionary -- Dictionary of secondary structure predictions for
            motifs and surrounding nucleotides per input file. The data is a 
            dictionary with the region as the key and the list of secondary 
            structure predictions as the value, arranged in three lists.
            {file_name: {"region_name": [structure, y_structure, gc_structure]}}. 
        output -- Output directory.
        secondary_per_region_check -- Check file for CSV file of secondary 
            structure data for each region.
            "Region_Name,Surrounding_Paired_Percent,Y_Paired_Percent,GC_Paired_Percent".
        secondary_per_file_check -- Check file for the CSV file of secondary 
            structure data for each sample.
            "File_Name,Surrounding_Paired_Percent,Y_Paired_Percent,GC_Paired_Percent".
        count_match_1 -- Integer, number of lines in the region check file.
        count_match_2 -- Integer, number of lines in the file check file.
     
"""

@pytest.mark.parametrize("file_dictionary, output, \
                          secondary_per_region_check, secondary_per_file_check, \
                          count_match_1, count_match_2", [
                        # Test 1: Multiple region in a single file.
                        # file_dictionary
                        ({"test1_file": {"name1": [['(','(','(','(','(','(','(','(','(','(','(','(','(','(','(','(','(','.'], 
                                                   ['(','(','(','(','.','.'], ['(','(','(','(','.','.']],
                                         "name2": [['.','.','.',')',')',')',')','(','(','(','(','(','(','(','(','(','(','(','(','(','(',')',')',')',')'], 
                                                   ['(',')','.','.',')',')','(',')'], ['(',')','.','.',')',')','.','.']]}},
                        # output
                        ("test_files/motif_secondary_features_structure/"),
                        # secondary_per_region_check
                        ("test_files/motif_secondary_features_structure/test1_secondary_per_region_check.csv"),
                        # secondary_per_file_check
                        ("test_files/motif_secondary_features_structure/test1_secondary_per_file_check.csv"),
                        # count_match_1
                        (3),
                        # count_match_2
                        (2)),
                        # Test 2: Multiple regions in multiple files.
                        # file_dictionary
                        ({"test1_file": {"name1": [['(','(','(','(','(','(','(','(','(','(','(','(','(','(','(','(','(','.'], 
                                                   ['(','(','(','(','.','.'], ['(','(','(','(','.','.']],
                                         "name2": [['.','.','.',')',')',')',')','(','(','(','(','(','(','(','(','(','(','(','(','(','(',')',')',')',')'], 
                                                   ['(',')','.','.',')',')','(',')'], ['(',')','.','.',')',')','.','.']]},
                          "test2_file": {"name1": [['(','(','(','(','(','(','(','(','(','(','(','(','(','(','(','(','(','.'], 
                                                   ['(','(','(','(','.','.'], ['(','(','(','(','.','.']],
                                         "name3": [['.','.','.','.',')',')',')','(','(','(','(','(','(','(','(','(','(','(','(','(','(',')',')',')',')'], 
                                                   ['(','.','.','.',')',')','(',')'], ['.',')','.','.',')',')','.','.']]},                         
                                                   },
                        # output
                        ("test_files/motif_secondary_features_structure/"),
                        # secondary_per_region_check
                        ("test_files/motif_secondary_features_structure/test2_secondary_per_region_check.csv"),
                        # secondary_per_file_check
                        ("test_files/motif_secondary_features_structure/test2_secondary_per_file_check.csv"),
                        # count_match_1
                        (4),
                        # count_match_2
                        (3)),            
                        ])

def test_output_region_and_file_data(file_dictionary, output,
                                 secondary_per_region_check, secondary_per_file_check,
                                 count_match_1, count_match_2):
    """
        GIVEN a file_dictionary with secondary structure predictions for
            motifs and surrounding nucleotides per input file.
        WHEN the function takes the input and outputs the secondary structure
            data for each region and each file.
        THEN the output files are checked for the correct data and titles.
    """
    # Executes function.
    m_s_f_s.output_region_and_file_data(file_dictionary, output)
    # Counts used to check output.
    line_count_1 = 0
    line_count_2 = 0
    # Checks to see if output directory ends with a "/".
    if output[-1] != "/":
        output = output + "/" 
    # Creates output file names.
    output_file_name_1 = (output + "secondary_structure_per_region.csv")
    output_file_name_2 = (output + "secondary_structure_per_file.csv")
    # Creates output files.
    output_file_1 = open(output_file_name_1, "r")
    output_file_2 = open(output_file_name_2, "r")
    # Checks each line of the output against the expected output
    # in the results_list. Counts should match.
    for line in output_file_1:
        # Cleans line to make it is easier to compare.
        line_clean = line.strip("\n")
        check1_opened = open(secondary_per_region_check, 'r')
        for line_check in check1_opened:
            # Cleans line to make it is easier to compare.
            line_check_clean = line_check.strip("\n")
            if line_clean == line_check_clean:
                line_count_1 += 1
    # Checks line counts.
    assert count_match_1 == line_count_1
    output_file_1.close()
    check1_opened.close()
    # Checks each line of the output against the expected output
    # in the results_list. Counts should match.
    for line in output_file_2:
        # Cleans line to make it is easier to compare.
        line_clean = line.strip("\n")
        check2_opened = open(secondary_per_file_check, 'r')
        for line_check in check2_opened:
            # Cleans line to make it is easier to compare.
            line_check_clean = line_check.strip("\n")
            if line_clean == line_check_clean:
                line_count_2 += 1
    # Checks line counts.
    assert count_match_2 == line_count_2
    output_file_2.close()
    check1_opened.close()   