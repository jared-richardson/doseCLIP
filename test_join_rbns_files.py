"""test_join_rbns_files.py- Unit tests for
    join_rbns_files.py
"""
import pytest
# Imports all functions from join_rbns_files.py.
import join_rbns_files as j_r_f

"""
Tests get_file_name with various file inputs.

    Data Types:
    file -- Input file string.
    file_name -- File name without 
        directory structure and 
        file extension.
"""

@pytest.mark.parametrize("file, file_name", [
                        # Test 1: File with full directory structure.
                        # file
                        (("/Users/username/test_files/random_enrichment_motif/test1.csv"),
                        # file_name
                        ("test1")),
                        # Test 2: File with no directory structure.
                        # file
                        (("test2.csv"),
                        # file_name
                        ("test2")),
                        ])

def test_get_file_name(file, file_name):
    """
        GIVEN a file input.
        WHEN the function parses the input file and returns the file name.
        THEN the output file name is checked for the correct name.
    """
    # Executes function and checks output.
    assert j_r_f.get_file_name(file) == file_name

"""
Tests join_and_output_rbns with various file inputs.

    Data Types:
    rbns_list -- (-r) List of CSV file of RBNS style motif enrichment values.
    output_file -- (-o) Output directory.   
    check_file -- CSV file of joined RBNS style motif enrichment values for 
        input files. First line contains file names used
        to generate motif enrichment values.
        File name- output_file + "all_joined_rbns_motif_enrichment.csv".
        Format - "Motif, Sample_Set1, Sample_Set2, ...".
    count_match -- Expected count of matching lines in output 
        and premade testing files.        
"""

@pytest.mark.parametrize("rbns_list, output_file, check_file, count_match", [
                        # Test 1: Three input files with identical motifs.
                        # rbns_list
                        ((["test_files/join_rbns_files/test_rbns_motif_enrichment.csv",
                           "test_files/join_rbns_files/test1_rbns_motif_enrichment.csv",
                           "test_files/join_rbns_files/test2_rbns_motif_enrichment.csv"]),
                        # output_file
                        ("test_files/join_rbns_files/test1"),
                        # check_file
                        ("test_files/join_rbns_files/test1/all_joined_rbns_motif_enrichment_check.csv"),
                        # count_match
                        (16)),
                        # Test 2: Three input files with different motifs and not all
                        # motifs present in all files.
                        # rbns_list
                        ((["test_files/join_rbns_files/test_rbns_motif_enrichment.csv",
                           "test_files/join_rbns_files/test3_rbns_motif_enrichment.csv",
                           "test_files/join_rbns_files/test4_rbns_motif_enrichment.csv"]),
                        # output_file
                        ("test_files/join_rbns_files/test2"),
                        # check_file
                        ("test_files/join_rbns_files/test2/all_joined_rbns_motif_enrichment_check.csv"),
                        # count_match
                        (18)),                     
                        ])

def test_join_and_output_rbns(rbns_list, output_file, check_file, count_match):
    """
        GIVEN a list of CSV files of RBNS style motif enrichment values that
            need to be joined and output in a directory.
        WHEN the function outputs a joned CSV file using the data from
            all the files in the required directory.
        THEN the output file is checked against a premade file to ensure
            the output is correct and the correct output location.
    """
    # Executes function.
    j_r_f.join_and_output_rbns(rbns_list, output_file)
    # Counts used to check output.
    line_count = 0
    # Checks to see if output directory ends with a "/".
    if output_file[-1] != "/":
        output_file = output_file + "/" 
    # Creates output file names.
    output_file_name = (output_file + "all_joined_rbns_motif_enrichment.csv")
    # Creates output files.
    output_file = open(output_file_name, "r")
    # Checks each line of the output against the expected output
    # in the results_list. Counts should match.
    for line in output_file:
        # Cleans line to make it is easier to compare.
        line_clean = line.strip("\n")
        check1_opened = open(check_file, 'r')
        for line_check in check1_opened:
            # Cleans line to make it is easier to compare.
            line_check_clean = line_check.strip("\n")
            if line_clean == line_check_clean:
                line_count += 1
    # Checks line counts.
    assert count_match == line_count
    output_file.close()
    check1_opened.close()