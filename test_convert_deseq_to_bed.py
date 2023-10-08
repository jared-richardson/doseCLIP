"""test_convert_deseq_to_bed.py- Unit tests for
    convert_deseq_to_bed.py
"""
import pytest
# Imports all functions from convert_deseq_to_bed.py.
import convert_deseq_to_bed as c_d_b

"""
Tests convert_regions() with different BED files.

Data Types:
    deseq_files -- List of DESeq2 produced files that need to be converted to BED format.
        Any DESeq2 file that contains gene names with the "~" separator should
        be able to be converted. Files will contain the same prefix as was input,
        and need to be in ".csv", ".txt", or ".tsv" format.
    bed_list -- List of BED files of the DESeq2 input are output.
    results_list -- List of results expected in function output. 
    count_match -- List of expected count of matching lines in output 
        and premade testing files.    
    """

@pytest.mark.parametrize("deseq_files, bed_list, \
                         results_list, count_match", [
                         # Test 1: One line DESeq2 filein CSV format.
                         # deseq_files
                         ((["test_files/convert_deseq_to_bed/deseq1.csv"]),
                          # bed_list
                          (["test_files/convert_deseq_to_bed/deseq1.bed"]),
                          # results_list
                          (["chr1_100_300_X_100_+"]),
                          # count_match  
                            ([1])),
                         # Test 2: Multiple line DESeq2 file in CSV format.
                         # deseq_files
                         ((["test_files/convert_deseq_to_bed/DESeq2.csv"]),
                          # bed_list
                          (["test_files/convert_deseq_to_bed/DESeq2.bed"]),
                          # results_list
                          (["chr1_100_300_X_100_+", "chr1_3100_3300_X_100_-",
                            "chr2_9000_9100_X_100_+", "chr4_1000_1100_X_100_+"]),
                          # count_match  
                            ([4])),
                         # Test 3: Multiple line DESeq2 file in TSV format with multiple
                         # file input.
                         # deseq_files
                         ((["test_files/convert_deseq_to_bed/deseq1.csv", 
                            "test_files/convert_deseq_to_bed/DESeq2.csv",
                            "test_files/convert_deseq_to_bed/deseq3.tsv"]),
                          # bed_list
                          (["test_files/convert_deseq_to_bed/deseq1.bed", 
                            "test_files/convert_deseq_to_bed/DESeq2.bed",
                            "test_files/convert_deseq_to_bed/deseq3.bed"]),
                          # results_list
                          (["chr1_100_300_X_100_+",
                            "chr1_100_300_X_100_+", "chr1_3100_3300_X_100_-",
                            "chr2_9000_9100_X_100_+", "chr4_1000_1100_X_100_+",
                            "chr7_100_300_X_100_+","chr7_100_300_X_100_-" ]),
                          # count_match  
                            ([1, 4, 2])),  
                         ])

def test_convert_regions(deseq_files, bed_list,
                         results_list, count_match):
    """
        GIVEN a DESeq2 file.
        WHEN the function parses and outputs each DESeq2 event
            and outputs it in BED format.
        THEN the output file is checked for the correct output by
            counting the number of perfectly matching lines.
    """
     # File count used to grab correct count_match
    file_count = 0
    # Executes function.
    c_d_b.convert_regions(deseq_files)
    # Opens testing output files to iterate through to check against
    # the expected output in the the results_list.
    for bed_file in bed_list:
        # File count used to grab correct count_match
        file_count 
        # Sets line_count to 0. This is the count of how many
        # lines match in the output and the expected output file.
        line_count = 0
        read1_opened = open(bed_file, 'r')
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