"""test_remove_duplicates_blast.py- Unit tests for
    remove_duplicates_blast.py.
"""
import pytest
# Imports all functions from remove_duplicates_blast.py.
import remove_duplicates_blast as r_d_b

"""Test remove_duplicates function using mutliple BLAST output files.

Date Types:
    blast_out -- List of BLAST produced output files.
        output -- Output directory.
    file_output_check -- List of CSV files of expected output.
    count_match_list  -- List of number of lines that match between the
        output file and the expected output file.
"""

@pytest.mark.parametrize("blast_out, output, file_output_check, count_match_list", [
                         # Test 1: Two line BLAST output file.
                         # blast_out
                         ((["test_files/remove_duplicates_blast/test_blast.tsv"]),
                         # output
                         ("test_files/remove_duplicates_blast"),
                         # file_output_check
                         (["test_files/remove_duplicates_blast/test_blast_single_region_check.tsv"]),
                         # count_match_list
                         ([1])),
                        # Test 2: Two BLAST files.
                        # blast_out
                        ((["test_files/remove_duplicates_blast/test_blast1.tsv",
                           "test_files/remove_duplicates_blast/test_blast2.tsv"]),
                        # output
                        ("test_files/remove_duplicates_blast"),
                        # file_output_check
                        (["test_files/remove_duplicates_blast/test_blast1_single_region_check.tsv",
                          "test_files/remove_duplicates_blast/test_blast2_single_region_check.tsv"]),
                          # count_match_list
                          ([2, 3]))    
                        ])

def test_remove_duplicates(blast_out, output, file_output_check, count_match_list):
    """
    GIVEN a list of BLAST output files that need duplicates removed.
    WHEN the function removces duplicates and outputs in the input
        output directory.
    THEN the output is checked for the correct filename and output.
    """
    # Runs remove_duplicates function.
    r_d_b.remove_duplicates(blast_out, output)
    # Checks for "/" in output directory.
    if output[-1] != "/":
        output += "/"
    # Loops through and checks output.
    for file in blast_out:
        # Removes ".tsv" from file name.    
        if file.find(".tsv") != -1:
            file = file.replace(".tsv", "")
        # Creates output file name.
        output_file = (output + file.split("/")[-1].split(".")[0] + "_single_region.tsv")
        # Count to keep track of index in count_match list.
        index_count = 0
        # Opens output file.
        with open(output_file, "r") as output_file:
            # Count used for matching lines.
            line_count = 0
            # Iterates through each line in the output file.
            for line in output_file:
                # Opens output check file.
                with open(file_output_check[index_count], "r") as file_output_check_file:
                    # Iterates through each line in the output check file.
                    for line_check in file_output_check_file:
                        # Checks if the lines match.
                        if line == line_check:
                            line_count += 1
                            pass
            # Checks if the number of lines that match is correct.
            assert line_count == count_match_list[index_count]
        # Adds one to index count.
        index_count += 1


