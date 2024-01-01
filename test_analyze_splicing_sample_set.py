"""test_analyze_splicing_sample_set.py- Unit tests for
    analyze_splicing_sample_set.py
"""
import pytest
# Imports all functions from analyze_splicing_sample_set.py.
import analyze_splicing_sample_set as a_s_s_s

"""Tests find_output_all_splicing() with three samples
    and different matching splicing events.

Data Types:
    filtered_splicing -- rMATS style output file that is
            filtered by delta PSI and an FDR of 10% or less.
            Contains the 'filtered_SE.MATS.JCEC.txt' suffix.
    output -- Prefix for output file.
    results_check -- Text file with results expected in function 
        produced files.    
    count_match -- Expected count of matching lines in output 
        and premade testing files.
    """

@pytest.mark.parametrize("filtered_splicing, output, results_check, \
                          count_match", [
                         # Test 1: Three common splicing events in
                         # three files. Some lines are only present
                         # in one/two files and identification fields
                         # should all be assessed.
                         # filtered_splicing
                         (('test_files/analyze_splicing_sample_set/test2_SE.MATS.JCEC.txt',
                           'test_files/analyze_splicing_sample_set/test3_SE.MATS.JCEC.txt',
                           'test_files/analyze_splicing_sample_set/test4_SE.MATS.JCEC.txt'),
                          # output
                          ("test_files/analyze_splicing_sample_set/test"),
                          # results_check
                          ("test_files/analyze_splicing_sample_set/test_all_filtered_SE.MATS.JCEC_check.txt"),
                          # count_match
                          (10)),
                         ])

def test_find_output_all_splicing(filtered_splicing, output, results_check, \
                                  count_match):
    """
        GIVEN a list of filtered splicing files that needs to be analyzed for
            overlapping splicing events in all files
        WHEN the function finds all common splicing events and outputs them
            with the filename added to the end of the splicing line.
        THEN the output file is compared against a premade file with the
            expected results and determined to be correct.
    """
    # Counts number of lines in the output file.
    line_count = 0
    a_s_s_s.find_output_all_splicing(filtered_splicing, output)
    # Opens output file for checking output.
    out1_opened = open(output + "_all_filtered_SE.MATS.JCEC.txt", "r")
    # Checks each line of the output against the expected output
    # in the results_list. Counts should match.
    for line in out1_opened:
        # Cleans line to make it is easier to compare.
        line_clean = line.strip("\n")
        check1_opened = open(results_check, 'r')
        for line_check in check1_opened:
            # Cleans line to make it is easier to compare.
            line_check_clean = line_check.strip("\n")
            if line_clean == line_check_clean:
                line_count += 1
    out1_opened.close()
    check1_opened.close()        
    # Checks line counts.
    assert count_match == line_count 

"""Tests find_output_multiple_binding() with various
    multiple binding regions per splicing events.

Data Types:
    splicing_binding -- CSV file with splicing and binding
        information for each region. Produced by the 
        match_binding_to_splicing.py script.
    output -- Prefix for output file.
    results_check -- Text file with results expected in function 
        produced files.    
    count_match_ -- Expected count of matching lines in output 
        and premade testing files.
"""

@pytest.mark.parametrize("splicing_binding, output, results_check, \
                         count_match", [
                        # Test 1: Single double region to splicing event
                        # with two files.
                        # splicing_binding
                        ((['test_files/analyze_splicing_sample_set/test1_splicing_binding.csv',
                          'test_files/analyze_splicing_sample_set/test12_splicing_binding.csv']),
                         # output
                         ("test_files/analyze_splicing_sample_set/test"),
                         # results_check_list
                         ("test_files/analyze_splicing_sample_set/test_multiple_binding_splicing_check.csv"),
                         # count_match_list
                         (3)),
                        # Test 2: Multiple double regions to splicing event with three files. 
                        # Tests each data field with multiple combinations between the three files.
                        # splicing_binding
                        ((['test_files/analyze_splicing_sample_set/test2_splicing_binding.csv',
                          'test_files/analyze_splicing_sample_set/test22_splicing_binding.csv',
                          'test_files/analyze_splicing_sample_set/test23_splicing_binding.csv']),
                         # output
                         ("test_files/analyze_splicing_sample_set/test2"),
                         # results_check_list
                         ("test_files/analyze_splicing_sample_set/test2_multiple_binding_splicing_check.csv"),
                         # count_match_list
                         (10)),
                        ])

def test_find_output_multiple_binding(splicing_binding, output, results_check, \
                                      count_match):
    """
        GIVEN a list of CSV files that have binding region to splicing
            information that needs to be analyzed for multiple binding
            regions per splicing event.
        WHEN the function finds all splicing events with multiple binding
            regions and outputs them per file.
        THEN the output file is compared against a premade file with the
            expected results and determined to be correct. 
    """
    a_s_s_s.find_output_multiple_binding(splicing_binding, output)
    # Counts used to check output.
    line_count = 0
    # Opens output file for checking output.
    out1_opened = open(output + "_multiple_binding_splicing.csv", "r")
    print(out1_opened)
    # Checks each line of the output against the expected output
    # in the results_list. Counts should match.
    for line in out1_opened:
        # Cleans line to make it is easier to compare.
        line_clean = line.strip("\n")
        check1_opened = open(results_check, 'r')
        for line_check in check1_opened:
            # Cleans line to make it is easier to compare.
            line_check_clean = line_check.strip("\n")
            if line_clean == line_check_clean:
                line_count += 1
    out1_opened.close()
    check1_opened.close()        
    # Checks line counts.
    assert count_match == line_count
