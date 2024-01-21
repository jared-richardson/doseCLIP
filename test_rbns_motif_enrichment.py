"""test_rbns_motif_enrichment.py- Unit tests for
    rbns_motif_enrichment.py
"""
import pytest
# Imports all functions from rbns_motif_enrichment.py.
import rbns_motif_enrichment as r_m_e

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
                        (("/Users/username/test_files/random_enrichment_motif/test1.fastq"),
                        # file_name
                        ("test1")),
                        # Test 2: File with no directory structure.
                        # file
                        (("test2.fastq"),
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
    assert r_m_e.get_file_name(file) == file_name

"""
Tests motif_detector with various FASTQ inputs.

    Data Types:
    fastq_file -- FASTQ file to be analyzed for motifs.
    kmer_size -- Size of kmer to be analyzed.
    motif_dictionary -- Dictionary of motif counts for each file.
        {"motif": normalized_count(float)}
"""

@pytest.mark.parametrize("fastq_file, kmer_size, motif_dictionary", [
                        # Test 1: FASTQ file with two sequences.
                        # fastq_file
                        (("test_files/rbns_motif_enrichment/test1.fastq"),
                        # kmer_size
                        (7),
                        # motif_dictionary
                        ({'TTTTTTT': 0.07142857142857142, 
                          'TTTTTTA': 0.03571428571428571, 
                          'TTTTTAA': 0.03571428571428571, 
                          'TTTTAAA': 0.03571428571428571, 
                          'TTTAAAA': 0.03571428571428571, 
                          'TTAAAAA': 0.03571428571428571, 
                          'TAAAAAA': 0.03571428571428571, 
                          'AAAAAAA': 0.03571428571428571, 
                          'TTTTTTG': 0.03571428571428571, 
                          'TTTTTGG': 0.03571428571428571, 
                          'TTTTGGG': 0.03571428571428571, 
                          'TTTGGGG': 0.03571428571428571, 
                          'TTGGGGG': 0.03571428571428571, 
                          'TGGGGGG': 0.03571428571428571, 
                          'GGGGGGG': 0.03571428571428571})),
                        # Test 2: FASTQ file with three sequences
                        # and N and n.
                        # fastq_file
                        (("test_files/rbns_motif_enrichment/test2.fastq"),
                        # kmer_size
                        (7),
                        # motif_dictionary
                        ({'TTTTTTT': 0.05357142857142857, 
                          'TTTTTTA': 0.03571428571428571, 
                          'TTTTTAA': 0.03571428571428571, 
                          'TTTTAAA': 0.03571428571428571, 
                          'TTTAAAA': 0.03571428571428571, 
                          'TTAAAAA': 0.03571428571428571, 
                          'TAAAAAA': 0.03571428571428571, 
                          'AAAAAAA': 0.03571428571428571, 
                          'TTTTTTG': 0.03571428571428571, 
                          'TTTTTGG': 0.03571428571428571, 
                          'TTTTGGG': 0.03571428571428571, 
                          'TTTGGGG': 0.03571428571428571, 
                          'TTGGGGG': 0.03571428571428571, 
                          'TGGGGGG': 0.03571428571428571, 
                          'GGGGGGG': 0.017857142857142856})),
                        ])

def test_motif_detector(fastq_file, kmer_size, motif_dictionary):
    """
        GIVEN a FASTQ file input and a kmer size for motif analysis.
        WHEN the function gets normalized enrichment values for each
            motif in the file.
        THEN the output dictionary is checked for the correct values.
    """
    # Executes function and checks output.
    assert r_m_e.motif_detector(fastq_file, kmer_size) == motif_dictionary

"""
Tests normalize_motifs with various file inputs.

    Data Types:
    motif_dictionary -- Dictionary of motif counts for each file.
        {"motif": normalized_count(float)}.
    motif_dictionary_sm -- Dictionary of motif counts for each file.
        {"motif": normalized_count(float)}.
    normalized_dictionary -- Dictionary of SM normalized motif counts
        for each file. {"motif": normalized_count(float)}.
"""

@pytest.mark.parametrize("motif_dictionary, motif_dictionary_sm, \
                         normalized_dictionary", [
                        # Test 1: Two dictionaries with multiple
                        # all matching motifs.
                        # motif_dictionary
                        (({'TTTTTTT': 3.0, 
                          'TTTTTTA': 4.0,
                          'AAAAAAA': 10.0,}),
                        # motif_dictionary_sm
                        ({'TTTTTTT': 1.0, 
                          'TTTTTTA': 8.0,
                          'AAAAAAA': 2.0,}),
                        # normalized_dictionary
                        ({'TTTTTTT': 3.0, 
                          'TTTTTTA': 0.5,
                          'AAAAAAA': 5.0,})),
                        # Test 1: Two dictionaries with one
                        # motif not present in both.
                        # motif_dictionary
                        (({'TTTTTTT': 3.0, 
                          'TTTTTAA': 4.0,
                          'AAAAAAA': 10.0,}),
                        # motif_dictionary_sm
                        ({'TTTTTTT': 1.0, 
                          'TTTTAAA': 8.0,
                          'AAAAAAA': 2.0,}),
                        # normalized_dictionary
                        ({'TTTTTTT': 3.0,
                          'AAAAAAA': 5.0,})),  
                        ])

def test_normalize_motifs(motif_dictionary, motif_dictionary_sm, \
                       normalized_dictionary):
    """
        GIVEN a two dictionaries with normalized motif
            enrichment values.
        WHEN the function additional normalizes the counts
            by dividing the first dictionary's value by the second.
        THEN the output dictionary is checked for the correct values.
    """
    # Executes function and checks output.
    assert r_m_e.normalize_motifs(motif_dictionary, 
                                  motif_dictionary_sm) == normalized_dictionary

"""
Tests get_motif_enrichment with various file inputs.

    Data Types:
    fastq_list -- (-f) List of paired-end FASTQ files of doseCLIP
        samples. Should contain all replicates and read1 and 
        read2 should be in same order.
    fastq_sm_list -- (-fs) List of paired-end FASTQ files of SM samples 
        pertaining to each and all input doseCLIP files. Should
        be in same order as doseCLIP samples, including with read1
        and read2. 
    kmer_size -- (-k) Size of kmer to be analyzed. 
    final_dictionary -- Dictionary of average motif enrichment
        values using all input files. 
        {"motif": normalized_count(float)}.
    file_dictionary -- Dictionary of motif counts for each file.
        {"file_name": {"motif": normalized_count(float)}}.
"""

@pytest.mark.parametrize("fastq_list, fastq_sm_list, kmer_size, \
                         final_dictionary, file_dictionary", [
                        # Test 1: Two FASTQ files input with identical
                        # motifs.
                        # fastq_list
                        ((["test_files/rbns_motif_enrichment/test1.fastq",
                           "test_files/rbns_motif_enrichment/test2.fastq"]),
                        # fastq_sm_list
                        (["test_files/rbns_motif_enrichment/test1.fastq",
                           "test_files/rbns_motif_enrichment/test2.fastq"]),
                        # kmer_size
                        (7),
                        # final_dictionary
                        ({'TTTTTTT': 1.0, 'TTTTTTA': 1.0, 
                          'TTTTTAA': 1.0, 'TTTTAAA': 1.0, 
                          'TTTAAAA': 1.0, 'TTAAAAA': 1.0, 
                          'TAAAAAA': 1.0, 'AAAAAAA': 1.0, 
                          'TTTTTTG': 1.0, 'TTTTTGG': 1.0, 
                          'TTTTGGG': 1.0, 'TTTGGGG': 1.0, 
                          'TTGGGGG': 1.0, 'TGGGGGG': 1.0, 
                          'GGGGGGG': 1.0}),
                        # file_dictionary
                        ({'test1-test1': {'TTTTTTT': 1.0, 'TTTTTTA': 1.0, 
                                    'TTTTTAA': 1.0, 'TTTTAAA': 1.0, 
                                    'TTTAAAA': 1.0, 'TTAAAAA': 1.0, 
                                    'TAAAAAA': 1.0, 'AAAAAAA': 1.0, 
                                    'TTTTTTG': 1.0, 'TTTTTGG': 1.0, 
                                    'TTTTGGG': 1.0, 'TTTGGGG': 1.0, 
                                    'TTGGGGG': 1.0, 'TGGGGGG': 1.0, 
                                    'GGGGGGG': 1.0}, 
                          'test2-test2': {'TTTTTTT': 1.0, 'TTTTTTA': 1.0, 
                                    'TTTTTAA': 1.0, 'TTTTAAA': 1.0, 
                                    'TTTAAAA': 1.0, 'TTAAAAA': 1.0, 
                                    'TAAAAAA': 1.0, 'AAAAAAA': 1.0, 
                                    'TTTTTTG': 1.0, 'TTTTTGG': 1.0, 
                                    'TTTTGGG': 1.0, 'TTTGGGG': 1.0, 
                                    'TTGGGGG': 1.0, 'TGGGGGG': 1.0, 
                                    'GGGGGGG': 1.0}})),
                        # Test 2: Two FASTQ files input with varying
                        # numbers of motifs.
                        # fastq_list
                        ((["test_files/rbns_motif_enrichment/test3.fastq",
                           "test_files/rbns_motif_enrichment/test4.fastq"]),
                        # fastq_sm_list
                        (["test_files/rbns_motif_enrichment/test1.fastq",
                          "test_files/rbns_motif_enrichment/test2.fastq"]),
                        # kmer_size
                        (7),
                        # final_dictionary
                        ({'TTTTTTT': 0.6507936507936508, 'TTTTTTA': 0.6190476190476191, 
                          'TTTTTAA': 0.6190476190476191, 'TTTTAAA': 0.6190476190476191, 
                          'TTTAAAA': 0.6190476190476191, 'TTAAAAA': 0.6190476190476191, 
                          'TAAAAAA': 0.6190476190476191, 'AAAAAAA': 0.5714285714285714, 
                          'TTTTTTG': 0.9523809523809523, 'TTTTTGG': 0.9523809523809523, 
                          'TTTTGGG': 0.9523809523809523, 'TTTGGGG': 0.9523809523809523, 
                          'TTGGGGG': 0.9523809523809523, 'TGGGGGG': 0.9523809523809523, 
                          'GGGGGGG': 1.5238095238095237}),
                        # file_dictionary
                        ({'test3-test1': {'TTTTTTT': 0.8571428571428572, 'TTTTTTA': 0.5714285714285714, 
                                          'TTTTTAA': 0.5714285714285714, 'TTTTAAA': 0.5714285714285714, 
                                          'TTTAAAA': 0.5714285714285714, 'TTAAAAA': 0.5714285714285714, 
                                          'TAAAAAA': 0.5714285714285714, 'AAAAAAA': 0.5714285714285714, 
                                          'TTTTTTG': 0.5714285714285714, 'TTTTTGG': 0.5714285714285714, 
                                          'TTTTGGG': 0.5714285714285714, 'TTTGGGG': 0.5714285714285714, 
                                          'TTGGGGG': 0.5714285714285714, 'TGGGGGG': 0.5714285714285714, 
                                          'GGGGGGG': 1.7142857142857144}, 
                          'test4-test2': {'TTTTTTT': 0.4444444444444444, 'TTTTTTG': 1.3333333333333333, 
                                          'TTTTTGG': 1.3333333333333333, 'TTTTGGG': 1.3333333333333333, 
                                          'TTTGGGG': 1.3333333333333333, 'TTGGGGG': 1.3333333333333333, 
                                          'TGGGGGG': 1.3333333333333333, 'GGGGGGG': 1.3333333333333333, 
                                          'TTTTTTA': 0.6666666666666666, 'TTTTTAA': 0.6666666666666666, 
                                          'TTTTAAA': 0.6666666666666666, 'TTTAAAA': 0.6666666666666666, 
                                          'TTAAAAA': 0.6666666666666666, 'TAAAAAA': 0.6666666666666666}})),
                        # Test 3: Three FASTQ files input with varying
                        # numbers of motifs.
                        # fastq_list
                        ((["test_files/rbns_motif_enrichment/test3.fastq",
                           "test_files/rbns_motif_enrichment/test4.fastq",
                           "test_files/rbns_motif_enrichment/test1.fastq"]),
                        # fastq_sm_list
                        (["test_files/rbns_motif_enrichment/test1.fastq",
                          "test_files/rbns_motif_enrichment/test2.fastq",
                          "test_files/rbns_motif_enrichment/test1.fastq"]),
                        # kmer_size
                        (7),
                        # final_dictionary
                        ({'TTTTTTT': 0.7671957671957671, 'TTTTTTA': 0.746031746031746, 
                          'TTTTTAA': 0.746031746031746, 'TTTTAAA': 0.746031746031746, 
                          'TTTAAAA': 0.746031746031746, 'TTAAAAA': 0.746031746031746, 
                          'TAAAAAA': 0.746031746031746, 'AAAAAAA': 0.7857142857142857, 
                          'TTTTTTG': 0.9682539682539683, 'TTTTTGG': 0.9682539682539683, 
                          'TTTTGGG': 0.9682539682539683, 'TTTGGGG': 0.9682539682539683, 
                          'TTGGGGG': 0.9682539682539683, 'TGGGGGG': 0.9682539682539683, 
                          'GGGGGGG': 1.349206349206349}),
                        # file_dictionary
                        ({'test3-test1': {'TTTTTTT': 0.8571428571428572, 'TTTTTTA': 0.5714285714285714, 
                                          'TTTTTAA': 0.5714285714285714, 'TTTTAAA': 0.5714285714285714, 
                                          'TTTAAAA': 0.5714285714285714, 'TTAAAAA': 0.5714285714285714, 
                                          'TAAAAAA': 0.5714285714285714, 'AAAAAAA': 0.5714285714285714, 
                                          'TTTTTTG': 0.5714285714285714, 'TTTTTGG': 0.5714285714285714, 
                                          'TTTTGGG': 0.5714285714285714, 'TTTGGGG': 0.5714285714285714, 
                                          'TTGGGGG': 0.5714285714285714, 'TGGGGGG': 0.5714285714285714, 
                                          'GGGGGGG': 1.7142857142857144}, 
                          'test4-test2': {'TTTTTTT': 0.4444444444444444, 'TTTTTTG': 1.3333333333333333, 
                                          'TTTTTGG': 1.3333333333333333, 'TTTTGGG': 1.3333333333333333, 
                                          'TTTGGGG': 1.3333333333333333, 'TTGGGGG': 1.3333333333333333, 
                                          'TGGGGGG': 1.3333333333333333, 'GGGGGGG': 1.3333333333333333, 
                                          'TTTTTTA': 0.6666666666666666, 'TTTTTAA': 0.6666666666666666, 
                                          'TTTTAAA': 0.6666666666666666, 'TTTAAAA': 0.6666666666666666, 
                                          'TTAAAAA': 0.6666666666666666, 'TAAAAAA': 0.6666666666666666}, 
                           'test1-test1': {'TTTTTTT': 1.0, 'TTTTTTA': 1.0, 'TTTTTAA': 1.0, 'TTTTAAA': 1.0, 
                                           'TTTAAAA': 1.0, 'TTAAAAA': 1.0, 'TAAAAAA': 1.0, 'AAAAAAA': 1.0, 
                                           'TTTTTTG': 1.0, 'TTTTTGG': 1.0, 'TTTTGGG': 1.0, 'TTTGGGG': 1.0, 
                                           'TTGGGGG': 1.0, 'TGGGGGG': 1.0, 'GGGGGGG': 1.0}})),                                 
                        ])

def test_get_motif_enrichment(fastq_list, fastq_sm_list, kmer_size,
                              final_dictionary, file_dictionary):
    """
        GIVEN a FASTQ file input and a kmer size for motif analysis.
        WHEN the function gets normalized enrichment values for each
            motif for all files.
        THEN the output dictionaries are checked for the correct values.
    """
    # Executes function and checks output.
    assert r_m_e.get_motif_enrichment(fastq_list, fastq_sm_list, 
                                      kmer_size) == (final_dictionary, 
                                                     file_dictionary)

"""
Tests output_motifs with various file inputs.

    Data Types:
    final_dictionary -- Dictionary of average motif enrichment
        values using all input files. 
        {"motif": normalized_count(float)}.
    file_dictionary -- Dictionary of motif counts for each file.
        {"file_name": {"motif": normalized_count(float)}}.
    output_file-- (-o) Output directory.
    file_prefix -- (-p) Prefix for output files.    
    check_file -- Check file for CSV file of RBNS style motif 
        enrichment values for input files. First line contains file names used
        to generate motif enrichment values.
        File name- output_file + file_prefix + "_rbns_motif_enrichment.csv".
        Format - "Motif, Normalized RBNS Motif Enrichment Value".
    count_match -- Expected count of matching lines in output 
        and premade testing files.        
"""

@pytest.mark.parametrize("final_dictionary, file_dictionary, output_file, \
                         file_prefix, check_file, count_match", [
                        # Test 1: Multiple motif output.
                        # final_dictionary
                        (({'TTTTTTT': 0.7671957671957671, 'TTTTTTA': 0.746031746031746, 
                          'TTTTTAA': 0.746031746031746, 'TTTTAAA': 0.746031746031746, 
                          'TTTAAAA': 0.746031746031746, 'TTAAAAA': 0.746031746031746, 
                          'TAAAAAA': 0.746031746031746, 'AAAAAAA': 0.7857142857142857, 
                          'TTTTTTG': 0.9682539682539683, 'TTTTTGG': 0.9682539682539683, 
                          'TTTTGGG': 0.9682539682539683, 'TTTGGGG': 0.9682539682539683, 
                          'TTGGGGG': 0.9682539682539683, 'TGGGGGG': 0.9682539682539683, 
                          'GGGGGGG': 1.349206349206349}),
                        # file_dictionary
                        ({'test3-test1': {'TTTTTTT': 0.8571428571428572, 'TTTTTTA': 0.5714285714285714, 
                                          'TTTTTAA': 0.5714285714285714, 'TTTTAAA': 0.5714285714285714, 
                                          'TTTAAAA': 0.5714285714285714, 'TTAAAAA': 0.5714285714285714, 
                                          'TAAAAAA': 0.5714285714285714, 'AAAAAAA': 0.5714285714285714, 
                                          'TTTTTTG': 0.5714285714285714, 'TTTTTGG': 0.5714285714285714, 
                                          'TTTTGGG': 0.5714285714285714, 'TTTGGGG': 0.5714285714285714, 
                                          'TTGGGGG': 0.5714285714285714, 'TGGGGGG': 0.5714285714285714, 
                                          'GGGGGGG': 1.7142857142857144}, 
                          'test4-test2': {'TTTTTTT': 0.4444444444444444, 'TTTTTTG': 1.3333333333333333, 
                                          'TTTTTGG': 1.3333333333333333, 'TTTTGGG': 1.3333333333333333, 
                                          'TTTGGGG': 1.3333333333333333, 'TTGGGGG': 1.3333333333333333, 
                                          'TGGGGGG': 1.3333333333333333, 'GGGGGGG': 1.3333333333333333, 
                                          'TTTTTTA': 0.6666666666666666, 'TTTTTAA': 0.6666666666666666, 
                                          'TTTTAAA': 0.6666666666666666, 'TTTAAAA': 0.6666666666666666, 
                                          'TTAAAAA': 0.6666666666666666, 'TAAAAAA': 0.6666666666666666}, 
                          'test1-test1': {'TTTTTTT': 1.0, 'TTTTTTA': 1.0, 'TTTTTAA': 1.0, 'TTTTAAA': 1.0, 
                                             'TTTAAAA': 1.0, 'TTAAAAA': 1.0, 'TAAAAAA': 1.0, 'AAAAAAA': 1.0, 
                                             'TTTTTTG': 1.0, 'TTTTTGG': 1.0, 'TTTTGGG': 1.0, 'TTTGGGG': 1.0, 
                                             'TTGGGGG': 1.0, 'TGGGGGG': 1.0, 'GGGGGGG': 1.0}}),
                        # output_file
                        ("test_files/rbns_motif_enrichment"),
                        # file_prefix
                        ("test"),
                        # check_file
                        ("test_files/rbns_motif_enrichment/test_rbns_motif_enrichment_check.csv"),
                        # count_match
                        (17)),                     
                        ])

def test_output_motifs(final_dictionary, file_dictionary, output_file,
                       file_prefix, check_file, count_match):
    """
        GIVEN a normalized motif enrichment dictionary and a file dictionary.
        WHEN the function outputs a CSV file of RBNS style motif enrichment
            values using the input dictionaries
        THEN the output file is compared against check_file for the correct
            number of matching lines.
    """
    # Executes function.
    r_m_e.output_motifs(final_dictionary, file_dictionary, output_file,
                        file_prefix)
    # Counts used to check output.
    line_count = 0
    # Checks to see if output directory ends with a "/".
    if output_file[-1] != "/":
        output_file = output_file + "/" 
    # Creates output file names.
    output_file_name = (output_file + file_prefix + "_rbns_motif_enrichment.csv")
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