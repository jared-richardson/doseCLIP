"""test_add_additional_sequence.py- Unit tests for
    add_additional_sequence.py
"""
import pytest
# Imports all functions from add_additional_sequence.py.
import add_additional_sequence as a_a_s

"""
Tests add_sequence_output() with different BED files.

Data Types:
    bed_list -- List of BED files of binding regions.
    length -- Length to add to both ends of all regions.
    output -- Output file directory.
    results_list -- List of results expected in function output. 
    count_match -- List of expected count of matching lines in output 
        and premade testing files.    
    """

@pytest.mark.parametrize("bed_list, length, output, \
                          results_list, count_match", [
                         # Test 1: One line BED file.
                         # bed_list
                         ((["test_files/add_additional_sequence/deseq1.bed"]),
                          # length
                          (500),
                          # output
                          ("test_files/add_additional_sequence"),
                          # results_list
                          (["chr1_1_800_X_100_+"]),
                          # count_match  
                            ([1])),
                         # Test 2: Multiple line BED file.
                         # bed_list
                         ((["test_files/add_additional_sequence/deseq2.bed"]),
                          # length
                          (500),
                          # output
                          ("test_files/add_additional_sequence/"),
                          # results_list
                          (["chr1_1_800_X_100_+", "chr1_2600_3800_X_100_-",
                            "chr2_8500_9600_X_100_+", "chr4_500_1600_X_100_+"]),
                          # count_match  
                          ([4])),
                         # Test 3: Multiple BED file input with multiple
                         # file input.
                         # bed_list
                         ((["test_files/add_additional_sequence/deseq31.bed", 
                            "test_files/add_additional_sequence/deseq32.bed",
                            "test_files/add_additional_sequence/deseq33.bed"]),
                          # length
                          (500),
                          # output
                          ("test_files/add_additional_sequence"),  
                          # results_list
                          (["chr1_1_800_X_100_+",
                            "chr1_1_800_X_100_+", "chr1_2600_3800_X_100_-",
                            "chr2_8500_9600_X_100_+", "chr4_500_1600_X_100_+",
                            "chr7_1_800_X_100_+","chr7_1_800_X_100_-" ]),
                          # count_match  
                            ([1, 4, 2])),  
                         ])

def test_add_sequence_output(bed_list, length, output, \
                             results_list, count_match):
    """
        GIVEN a BED file that needs additional sequence added to the start
            and end of each region.
        WHEN the function adds additional sequence to the start and end of
            each region.
        THEN the output file is checked for the correct output by
            counting the number of perfectly matching lines.
    """
     # File count used to grab correct count_match.
    file_count = 0
    # Executes function.
    a_a_s.add_sequence_output(bed_list, length, output)
    # Opens testing output files to iterate through to check against
    # the expected output in the the results_list.
    for bed in bed_list:
        # File count used to grab correct count_match.
        file_count 
        # Sets line_count to 0. This is the count of how many
        # lines match in the output and the expected output file.
        line_count = 0
        if output[-1] == "/":
            output_file_name = (output + bed.split("/")[-1].split(".")[0] + \
                                "_additional_sequence.bed")
        else:
            output_file_name = (output + "/" + bed.split("/")[-1].split(".")[0] + \
                                "_additional_sequence.bed") 
        read1_opened = open(output_file_name, 'r')
        # Checks each line of the output against the expected output
        # in the results_list. Counts should match.
        for line in read1_opened:
            # Makes line so it is easier to compare.
            line_clean = line.strip("\n").replace("\t","_")
            if line_clean in results_list:
                line_count += 1
        # Checks line counts.
        assert count_match[file_count] == line_count
        file_count += 1
        read1_opened.close()