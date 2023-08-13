"""test_join_binding_regions.py- Unit tests for
    join_binding_regions.py
"""
import pytest
# Imports all functions from join_binding_regions.py.
import count_binding_regions as c_b_r

"""
Tests count_binding_regions() with different DeSeq2 input files,
    different BED files, and different overlapping events.

Data Types:
    deseq2_files -- List of DeSeq2 produced differential expression
        file(s). Must be CSV file.
    region_bed_file -- BED file with binding regions. This file
        would be produced from a separate experiment where there
        where no SM Input samples.
    counts_file_csv -- Output counts file. Contains overlapping
        event counts for each input file.
    results_list -- List of results expected in function produced joined_gtf_file. 
    count_match -- Expected count of matching lines in output and premade
        testing files.    
    """

@pytest.mark.parametrize("deseq2_files, region_bed_file, counts_file_csv, \
                         results_list, count_match", [
                         # Test 1: One overlapping DeSeq2 and BED event.
                         # deseq2_files
                         ((["test_files/count_binding_regions/deseq2_differential_expression1.csv"]),
                          # region_bed_file
                          ("test_files/count_binding_regions/test1.bed"),
                          # counts_file_csv
                          ("test_files/count_binding_regions/counts_file1.csv"),
                          # results_list
                          (["DeSeq2_filename,BED_filename,Overlapping_count",
                            "deseq2_differential_expression1.csv,test1.bed,1"]),
                          # count_match  
                            (2)),
                         # Test 2: Multiple overlapping DeSeq2 and BED events.
                         # deseq2_files
                         ((["test_files/count_binding_regions/deseq2_differential_expression2.csv"]),
                          # region_bed_file
                          ("test_files/count_binding_regions/test2.bed"),
                          # counts_file_csv
                          ("test_files/count_binding_regions/counts_file2.csv"),
                          # results_list
                          (["DeSeq2_filename,BED_filename,Overlapping_count",
                            "deseq2_differential_expression2.csv,test2.bed,7"]),
                          # count_match  
                            (2)),
                         # Test 3: Multiple overlapping DeSeq2 and BED events and
                         # multiple DeSeq2 input files.
                         # deseq2_files
                         ((["test_files/count_binding_regions/deseq2_differential_expression2.csv",
                            "test_files/count_binding_regions/deseq2_differential_expression3.csv",
                            "test_files/count_binding_regions/deseq2_differential_expression4.csv"]),
                          # region_bed_file
                          ("test_files/count_binding_regions/test2.bed"),
                          # counts_file_csv
                          ("test_files/count_binding_regions/counts_file2.csv"),
                          # results_list
                          (["DeSeq2_filename,BED_filename,Overlapping_count",
                            "deseq2_differential_expression2.csv,test2.bed,7",
                            "deseq2_differential_expression2.csv,test2.bed,7",
                            "deseq2_differential_expression2.csv,test2.bed,3"]),
                          # count_match  
                            (2)),     
                         ])

def test_count_binding_regions(deseq2_files, region_bed_file, counts_file_csv, \
                               results_list, count_match):
    """
        GIVEN a DeSeq2 produced file or files, and a Piranha produced BED file
            from another CLIP experiment that need to be compared for overlapping
            events.
        WHEN the function function finds overlapping binding events between
            all input files, counts the events, and outputs the data
        THEN the output file is checked for the correct output by
            counting the number of perfectly matching lines.
    """
    # Sets line_count to 0. This is the count of how many
    # lines match in the output and the expected output file.
    line_count = 0
    # Executes function.
    c_b_r.count_binding_regions(deseq2_files, region_bed_file, counts_file_csv)
    # Opens testing output files to iterate through to check against
    # the expected output in the the results_list.
    read1_opened = open(counts_file_csv, 'r')
    # Checks each line of the output against the expected output
    # in the results_list. Counts should match.
    for line in read1_opened:
        # Cleans line to make it is easier to compare.
        line_clean = line.strip("\n").replace("\t","_")
        if line_clean in results_list:
            line_count += 1
    # Checks line counts.
    assert count_match == line_count