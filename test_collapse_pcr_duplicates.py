"""test_collapse_pcr_duplicates.py- Unit tests for
    collapse_pcr_duplicates.py
"""
import pytest
# Imports all functions from collapse_pcr_duplicates.py.
import collapse_pcr_duplicates as c_p_d

"""
Tests reverse_complement() with different kmers.

Data Types:
    kmer- Takes a kmer of any length and returns the
        reverse complement of the nucleotides.
    reverse_complement- Reverse complement of
        input sequence. This is the output of the
        function.
    """

@pytest.mark.parametrize("kmer, reverse_complement", [
                         # Test 1: Normal kmer input.
                         # kmer
                         (("ATGC"),
                          # reverse_complement
                          ("GCAT")),
                         # Test 2: Long kmer input.
                         # kmer
                         (("ATGCGCCCTAGT"),
                          # reverse_complement
                          ("ACTAGGGCGCAT")),
                         # Test 3: Short kmer input.
                         # kmer
                         (("A"),
                          # reverse_complement
                          ("T")),
                         # Test 4: None Nucleotide kmer input.
                         # kmer
                         (("ATGCFR"),
                          # reverse_complement
                          ("RFGCAT")),
                         # Test 5: Empty input.
                         # kmer
                         ((""),
                          # reverse_complement
                          (""))
                         ])

def test_reverse_complement(kmer, reverse_complement):
    """
        GIVEN a kmer that needs to be reverse complemented.
        WHEN the function reverse complements the kmer.
        THEN the reverse complement is checked for correctness.
    """
    assert c_p_d.reverse_complement(kmer) == reverse_complement


"""
Tests deduplicate_fastq() with different kmers.

Data Types:
    fastq_read1- First read of paired-end sequencing data 
        (FASTQ format).
    fastq_read2- Second read of paired-end sequencing data 
        (FASTQ format).   
    fastq_read1_out- First read expected output FASTQ file 
        with PCR duplicates removed.
    fastq_read2_out- Second read expected output FASTQ file
        with PCR duplicates removed.
    count_match- Expected count of matching lines in output and premade
        testing files.
    """

@pytest.mark.parametrize("fastq_read1, fastq_read2, \
                          fastq_read1_out, fastq_read2_out, \
                          count_match", [
                         # Test 1: Normal FASTQ input (two identifcal reads).
                         # fastq_read1
                         (("test_files/collapse_pcr_duplicates/test1_read1.fastq"),
                          # fastq_read2
                          ("test_files/collapse_pcr_duplicates/test1_read2.fastq"),
                          # fastq_read1_out
                          ("test_files/collapse_pcr_duplicates/test1_read1_check.fastq"),
                          # fastq_read2_out
                          ("test_files/collapse_pcr_duplicates/test1_read2_check.fastq"),
                          # count_match
                          (4)),
                         # Test 2: Normal FASTQ input (three identical reads).
                         # fastq_read1
                         (("test_files/collapse_pcr_duplicates/test2_read1.fastq"),
                          # fastq_read2
                          ("test_files/collapse_pcr_duplicates/test2_read2.fastq"),
                          # fastq_read1_out
                          ("test_files/collapse_pcr_duplicates/test2_read1_check.fastq"),
                          # fastq_read2_out
                          ("test_files/collapse_pcr_duplicates/test2_read2_check.fastq"),
                          # count_match
                          (4)),
                         # Test 3: Normal FASTQ input (second read of one input does not match, three other identical reads).
                         # fastq_read1
                         (("test_files/collapse_pcr_duplicates/test3_read1.fastq"),
                          # fastq_read2
                          ("test_files/collapse_pcr_duplicates/test3_read2.fastq"),
                          # fastq_read1_out
                          ("test_files/collapse_pcr_duplicates/test3_read1_check.fastq"),
                          # fastq_read2_out
                          ("test_files/collapse_pcr_duplicates/test3_read2_check.fastq"),
                          # count_match
                          (8)),
                         # Test 4: Normal FASTQ input (two sets of matching reads, one read with no "adaptor1", and a second read UMI that
                         # does not match the first).
                         # fastq_read1
                         (("test_files/collapse_pcr_duplicates/test4_read1.fastq"),
                          # fastq_read2
                          ("test_files/collapse_pcr_duplicates/test4_read2.fastq"),
                          # fastq_read1_out
                          ("test_files/collapse_pcr_duplicates/test4_read1_check.fastq"),
                          # fastq_read2_out
                          ("test_files/collapse_pcr_duplicates/test4_read2_check.fastq"),
                          # count_match
                          (16))                             
                         ])

def test_deduplicate_fastq(tmpdir, fastq_read1, fastq_read2,
                            fastq_read1_out, fastq_read2_out, count_match):
    """
        GIVEN a set of FASTQ files that need PCR duplicates removed.
        WHEN the function removes PCR duplicates and removes UMIs.
        THEN the temporary output FASTQ is checked for the proper deduplicated
            reads and sequences against what is expected.
    """
    # Sets line_count to 0. This is the count of how many
    # lines match in the output and the expected output file.
    line_count = 0
    # Opens testing output files to iterate through to check against
    # temparory output below.
    read1_opened = open(fastq_read1_out, 'r')
    read2_opened = open(fastq_read2_out, 'r')
    # Makes temporary output files.
    fastq_temp_1 = tmpdir.join("fastq_temp_1")
    fastq_temp_2 = tmpdir.join("fastq_temp_2")
    # Executes function.
    c_p_d.deduplicate_fastq(fastq_read1, fastq_read2, fastq_temp_1.strpath, fastq_temp_2.strpath)
    # Checks each line of the temporary output against the expected output
    # in the premade test files. Counts should match.
    for line in read1_opened:
        if (fastq_temp_1.read().find(line) > -1):
            line_count += 1
    # Checks line counts.
    assert count_match == line_count
    # Resets line_count for read 2 file.
    line_count = 0
    # Performs test for read 2.
    for line in read2_opened:
        if (fastq_temp_2.read().find(line) > -1):
            line_count += 1
    # Checks line counts.
    assert count_match == line_count              
    