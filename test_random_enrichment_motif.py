"""test_random_enrichment_motif.py- Unit tests for
    random_enrichment_motif.py
"""
import pytest
# Imports all functions from make_random_bed_for_motif.py.
import random_enrichment_motif as r_e_m

"""
Tests get_average_random with various motif_dictionary 
    inputs and numbers of dictionaries.

    Data Types:
    motif_dictionary -- Dictionary of motif counts for each file.
        {"file_name": {"motif": normalized_count(float)}}
    random_dictionary -- Output Dictionary of average motif counts 
        for all random files
        {"motif": normalized_count(float)}
"""

@pytest.mark.parametrize("motif_dictionary, random_dictionary", [
                        # Test 1: Two random input dictionaries with single input.
                        # motif_dictionary
                        (({"file_random1": {"motif1": 1.0},
                           "file_random2": {"motif1": 1.0}}),
                        # random_dictionary
                        ({"motif1": 1.0})),
                        # Test 2: Two random input dictionaries with multiple inputs.
                        # motif_dictionary
                        (({"file_random1": {"motif1": 1.0, "motif2": 2.0},
                           "file_random2": {"motif1": 1.0, "motif2": 4.0}}),
                        # random_dictionary
                        ({"motif1": 1.0, "motif2": 3.0})),
                        # Test 3: Three random input dictionaries with multiple inputs.
                        # motif_dictionary
                        (({"file_random1": {"motif1": 1.0, "motif2": 2.0, "motif3": 3.0},
                           "file_random2": {"motif1": 1.0, "motif2": 4.0, "motif3": 6.0},
                           "file_random3": {"motif1": 1.0, "motif2": 6.0}}),
                        # random_dictionary
                        ({"motif1": 1.0, "motif2": 4.0, "motif3": 3.0})),
                        ])

def test_get_average_random(motif_dictionary, random_dictionary):
    """
        GIVEN a motif_dictionary of random motif counts for each 
            random control file.
        WHEN the function averages motif counts for all files and
            organizes the counts by motif in a dictionary.
        THEN the output dictionary is checked for the correct entries.
    """
    # Executes function.
    output_dictionary = r_e_m.get_average_random(motif_dictionary)
    # Checks output_dictionary against random_dictionary.
    assert output_dictionary == random_dictionary

"""
Tests motif_detector with various FASTA file inputs and motif
    sequence inputs.

    Data Types:
    region_fasta_list -- List of FASTA files of RBP binding regions.
            Should include previously generated random regions from
            make_random_bed_for_motif.py. The random files are detected
            by searching for "random" in the file name.
    kmer_size -- Size of kmer to be analyzed.
    output_file_prefix -- Prefix for output files.
    motif_list -- List of motifs to be analyzed. Default is YGCY.
    motif_title -- Title for motif output file. Default is YGCY.   
    file_dictionary -- Dictionary of motif counts for each file.
        {"file_name": {"motif": normalized_count(float)}}
    region_dictionary -- Dictionary of motif counts for each region.
        {"event_name": normalized_motif_count(float)}
    motif_title -- Title for motif output file. Default is YGCY. 
"""

@pytest.mark.parametrize("region_fasta_list, kmer_size, motif_list, motif_title, \
                         file_dictionary, region_dictionary", [
                        # Test 1: Single region FASTA file with default YGCY motif.
                        # Also uses default motif_title and motif_list.
                        # region_fasta_list
                        ((["test_files/random_enrichment_motif/test1.fasta"]),
                        # kmer_size
                        (4),
                        # motif_list
                        (),
                        # motif_title
                        (),
                        # file_dictionary
                        ({"test1.fasta": {'CCGC': 0.125,
                                          'CGCC': 0.125,
                                          'CGCT': 0.125,
                                          'GCCG': 0.125,
                                          'TCGC': 0.125}}),
                        # region_dictionary
                        ({"name1": 0.25})),
                        # Test 2: Multiple region FASTA file with default YGCY motif.
                        # Default motif_title and motif_list are manually input. Also
                        # checks for all YGCY variants.
                        # region_fasta_list
                        ((["test_files/random_enrichment_motif/test2.fasta"]),
                         # kmer_size
                         (4),
                         # motif_list
                         (["YGCY"]),
                         # motif_title
                         ("YGCY"),
                         # file_dictionary
                         ({"test2.fasta": {'TCGC': 0.1,
                                           'CGCC': 0.1,
                                           'GCCG': 0.05,
                                           'CCGC': 0.05,
                                           'CGCT': 0.1,
                                           'GCTG': 0.05,
                                           'CTGC': 0.05,
                                           'TGCT': 0.05}}),
                         # region_dictionary
                         ({"name1": 0.25,
                           "name2": 0.25,
                           "name3": 0.25})),
                        # Test 3: Multiple region FASTA file with custom input
                        # motifs.
                        # region_fasta_list
                        ((["test_files/random_enrichment_motif/test3.fasta"]),
                         # kmer_size
                         (4),
                         # motif_list
                         (["TTTT", "CCCC"]),
                         # motif_title
                         ("TTTTCCCC"),
                         # file_dictionary
                         ({"test3.fasta": {'TTTT': 0.07142857142857142,
                                           'TTTC': 0.07142857142857142,
                                           'CCCC': 0.07142857142857142,
                                           'CCCT': 0.07142857142857142,
                                           'CGCC': 0.07142857142857142}}),
                         # region_dictionary
                         ({"name1": 0.2,
                           "name2": 0.2,
                           "name3": 0.00})),
                        # Test 4: Multiple region FASTA file with default YGCY motif.
                        # Multiple input files
                        # region_fasta_list
                        ((["test_files/random_enrichment_motif/test1.fasta",
                           "test_files/random_enrichment_motif/test2.fasta",
                           "test_files/random_enrichment_motif/test4.fasta"]),
                         # kmer_size
                         (4),
                         # motif_list
                         (["YGCY"]),
                         # motif_title
                         ("YGCY"),
                         # file_dictionary
                         ({"test1.fasta": {'CCGC': 0.125,
                                          'CGCC': 0.125,
                                          'CGCT': 0.125,
                                          'GCCG': 0.125,
                                          'TCGC': 0.125},
                           "test2.fasta": {'TCGC': 0.1,
                                           'CGCC': 0.1,
                                           'GCCG': 0.05,
                                           'CCGC': 0.05,
                                           'CGCT': 0.1,
                                           'GCTG': 0.05,
                                           'CTGC': 0.05,
                                           'TGCT': 0.05},
                           "test4.fasta": {'TCGC': 0.1,
                                           'CGCC': 0.1,
                                           'GCCG': 0.05,
                                           'CCGC': 0.05,
                                           'CGCT': 0.1,
                                           'GCTG': 0.05,
                                           'CTGC': 0.05,
                                           'TGCT': 0.05}}),
                         # region_dictionary
                         ({"name1": 0.25,
                           "name2": 0.25,
                           "name3": 0.25,
                           "name41": 0.25,
                           "name42": 0.25,
                           "name43": 0.25})),   
                        ])

def test_motif_detector(region_fasta_list, kmer_size, motif_list, motif_title,
                        file_dictionary, region_dictionary):
    """
        GIVEN a list of FASTA files of RBP binding regions.
        WHEN the function analyzes each region for motifs and
            outputs the motif counts for each region and file.
        THEN the output dictionaries are checked for the correct entries.
    """
    # Executes function.
    if len(motif_list) == 0:
        (output_file_dictionary, output_region_dictionary, 
         motif_list, motif_title) = r_e_m.motif_detector(region_fasta_list, kmer_size)
    else:    
        (output_file_dictionary, output_region_dictionary,
         motif_list, motif_title) = r_e_m.motif_detector(region_fasta_list, kmer_size,
                                                         motif_list, motif_title)
    # Checks output_motif_dictionary against motif_dictionary.
    assert output_file_dictionary == file_dictionary
    # Checks output_region_dictionary against region_dictionary.
    assert output_region_dictionary == region_dictionary

"""
Tests output_motifs with various FASTA file inputs and motif
    sequence inputs.

    Data Types:
    file_dictionary -- Dictionary of motif counts for each file.
        {"file_name": {"motif": normalized_count(float)}}
    region_dictionary -- Dictionary of motif counts for each region.
        {"event_name": normalized_motif_count(float)}
    output_file_prefix -- Prefix for output files.
    motif_list -- List of motifs to be analyzed. Default is YGCY.
    motif_title -- Title for motif output file. Default is YGCY.
    output_region_csv -- CSV files of random motif enrichment values
        for each binding region from all input files. File name-
        output_file_prefix + "_motifs_per_region.csv". Format -
        Region Name, Random Motif Enrichment Value.
    output motif CSV -- CSV file of random motif enrichment values
        per target motif per file. Last section indicates whether
        motif was able to be normalized using the random files.
        All should be normalized but was added for edge cases.
        File name- output_file_prefix + "total_motifs.csv". 
        Format -
        File Name, Motif, Random Motif Enrichment Value, Normalized.
    output all motifs CSV -- CSV file of random motif enrichment values
        for all motifs per file. Last section indicates whether
        motif was able to be normalized using the random files.
        All should be normalized but was added for edge cases.
        File name- output_file_prefix + "all_motifs.csv". 
        Format -
        File Name, Motif, Random Motif Enrichment Value, Normalized.
    results_dictionary -- Dictionary of expected results from each 
        output file.
    count_match -- List of expected count of matching lines in output 
        and premade testing files.  
"""

@pytest.mark.parametrize("file_dictionary, region_dictionary, output_file_prefix, \
                          motif_list, motif_title, \
                          output_region_csv, output_motif_csv, \
                          output_all_motifs_csv, results_dictionary, count_match", [
                        # Test 1: Single region FASTA file with default YGCY motif.
                        # Also uses default motif_title and motif_list. One identical
                        # random file is used.
                        # file_dictionary
                        (({"test1.fasta": {'CCGC': 0.125,
                                          'CGCC': 0.125,
                                          'CGCT': 0.125,
                                          'GCCG': 0.125,
                                          'TCGC': 0.125},
                          "test1_random.fasta": {'CCGC': 0.125,
                                                  'CGCC': 0.125,
                                                  'CGCT': 0.125,
                                                  'GCCG': 0.125,
                                                  'TCGC': 0.125}}),
                        # region_dictionary
                        ({"name1": 0.25}),
                        # output_file_prefix
                        ("test_files/random_enrichment_motif/test1"),
                        # motif_list
                        (["CGCC", "CGCT", "TGCC", "TGCT"]),
                        # motif_title
                        ("YGCY"),                 
                        # output_region_csv
                        ("test_files/random_enrichment_motif/test1_motifs_per_region.csv"),
                        # output_motif_csv
                        ("test_files/random_enrichment_motif/test1_total_motifs.csv"),
                        # output_all_motifs_csv
                        ("test_files/random_enrichment_motif/test1_all_motifs.csv"),
                        # results_dictionary
                        ({"test_files/random_enrichment_motif/test1_motifs_per_region.csv": \
                            ["Region Name, YGCY Motif Per 100 Nucleotides Value",
                             "name1,25.0"],
                          "test_files/random_enrichment_motif/test1_total_motifs.csv": \
                            ["File Name, YGCY Motif, Random Motif Enrichment Value, Normalized",
                             "test1.fasta,CGCC,1.0,Yes",
                             "test1.fasta,CGCT,1.0,Yes"],
                           "test_files/random_enrichment_motif/test1_all_motifs.csv": \
                            ["File Name, Motif, Random Motif Enrichment Value, Normalized",
                             "test1.fasta,CCGC,1.0,Yes",
                             "test1.fasta,CGCC,1.0,Yes",
                             "test1.fasta,CGCT,1.0,Yes",
                             "test1.fasta,GCCG,1.0,Yes",
                             "test1.fasta,TCGC,1.0,Yes"]}),
                        # count_match
                        ([2, 3, 6])),
                        # Test 2: Multiple FASTA file input with multuple random files.
                        # file_dictionary
                        (({"test1.fasta": {'CCGC': 0.125,
                                          'CGCC': 0.125,
                                          'CGCT': 0.125,
                                          'GCCG': 0.125,
                                          'TCGC': 0.125},
                           "test2.fasta": {'TCGC': 0.1,
                                           'CGCC': 0.1,
                                           'GCCG': 0.05,
                                           'CCGC': 0.05,
                                           'CGCT': 0.1,
                                           'GCTG': 0.05,
                                           'CTGC': 0.05,
                                           'TGCT': 0.05},
                           "test4.fasta": {'TCGC': 0.1,
                                           'CGCC': 0.1,
                                           'GCCG': 0.05,
                                           'CCGC': 0.05,
                                           'CGCT': 0.1,
                                           'GCTG': 0.05,
                                           'CTGC': 0.05,
                                           'TGCT': 0.0},                
                          "test1_random.fasta": {'CCGC': 0.125,
                                                  'CGCC': 0.125,
                                                  'CGCT': 0.125,
                                                  'GCCG': 0.125,
                                                  'TCGC': 0.125},
                          "test2_random.fasta": {'CCGC': 0.375,
                                                  'CGCC': 0.375,
                                                  'CGCT': 0.375,
                                                  'GCCG': 0.375,
                                                  'TCGC': 0.375},
                          "test3_random.fasta": {'CCGC': 0.5}}),
                        # region_dictionary
                        ({"name1": 0.25,
                           "name2": 0.25,
                           "name3": 0.25,
                           "name41": 0.25,
                           "name42": 0.50,
                           "name43": 0.50}),
                        # output_file_prefix
                        ("test_files/random_enrichment_motif/test2"),
                        # motif_list
                        (["CGCC", "CGCT", "TGCC", "TGCT"]),
                        # motif_title
                        ("YGCY"),                 
                        # output_region_csv
                        ("test_files/random_enrichment_motif/test2_motifs_per_region.csv"),
                        # output_motif_csv
                        ("test_files/random_enrichment_motif/test2_total_motifs.csv"),
                        # output_all_motifs_csv
                        ("test_files/random_enrichment_motif/test2_all_motifs.csv"),
                        # results_dictionary
                        ({"test_files/random_enrichment_motif/test2_motifs_per_region.csv": \
                            ["Region Name, YGCY Motif Per 100 Nucleotides Value",
                             "name1,25.0",
                             "name2,25.0",
                             "name3,25.0",
                             "name41,25.0",
                             "name42,50.0",
                             "name43,50.0"],
                          "test_files/random_enrichment_motif/test2_total_motifs.csv": \
                            ["File Name, YGCY Motif, Random Motif Enrichment Value, Normalized",
                             "test1.fasta,CGCC,0.75,Yes",
                             "test1.fasta,CGCT,0.75,Yes",
                             "test2.fasta,CGCC,0.6000000000000001,Yes",
                             "test2.fasta,CGCT,0.6000000000000001,Yes",
                             "test2.fasta,TGCT,0.05,No",
                             "test4.fasta,CGCC,0.6000000000000001,Yes",
                             "test4.fasta,CGCT,0.6000000000000001,Yes",
                             "test4.fasta,TGCT,0.0,No"],
                           "test_files/random_enrichment_motif/test2_all_motifs.csv": \
                            ["File Name, Motif, Random Motif Enrichment Value, Normalized",
                             "test1.fasta,CCGC,0.375,Yes",
                             "test1.fasta,CGCC,0.75,Yes",
                             "test1.fasta,CGCT,0.75,Yes",
                             "test1.fasta,GCCG,0.75,Yes",
                             "test1.fasta,TCGC,0.75,Yes",
                             "test2.fasta,TCGC,0.6000000000000001,Yes",
                             "test2.fasta,CGCC,0.6000000000000001,Yes",
                             "test2.fasta,GCCG,0.30000000000000004,Yes",
                             "test2.fasta,CCGC,0.15000000000000002,Yes",
                             "test2.fasta,CGCT,0.6000000000000001,Yes",
                             "test2.fasta,GCTG,0.05,No",
                             "test2.fasta,CTGC,0.05,No",
                             "test2.fasta,TGCT,0.05,No",
                             "test4.fasta,TCGC,0.6000000000000001,Yes",
                             "test4.fasta,CGCC,0.6000000000000001,Yes",
                             "test4.fasta,GCCG,0.30000000000000004,Yes",
                             "test4.fasta,CCGC,0.15000000000000002,Yes",
                             "test4.fasta,CGCT,0.6000000000000001,Yes",
                             "test4.fasta,GCTG,0.05,No",
                             "test4.fasta,CTGC,0.05,No",
                             "test4.fasta,TGCT,0.0,No"]}),
                        # count_match
                        ([7, 9, 22])),
                       ]) 

def test_output_motifs(file_dictionary, region_dictionary, output_file_prefix, motif_list,
                          motif_title, output_region_csv, output_motif_csv,
                          output_all_motifs_csv, results_dictionary, count_match):
    """
        GIVEN a dictionary of motif counts for each file and region.
        WHEN the function outputs the motif counts for each region and file.
        THEN the output files are checked for the correct entries.
    """
    # Executes function.
    r_e_m.output_motifs(file_dictionary, region_dictionary,
                        output_file_prefix, motif_list,
                        motif_title)
    # Names and opens output files from function.
    output_region_csv_open = open(output_region_csv, 'r')
    output_motif_csv_open = open(output_motif_csv, 'r')
    output_all_motifs_csv_open = open(output_all_motifs_csv, 'r')
    # List used to save matching counts.
    count_list = []
    # Count used to count output that matches expected results.
    line_count = 0
    # Checks output in each file and compares to expected results.
    region_list = results_dictionary[output_region_csv]
    # Iterates through each matching file and checks for
    # correct output.
    for line in output_region_csv_open:
        line_clean = line.strip("\n")
        if line_clean in region_list:
            line_count += 1
    count_list.append(line_count)
    output_motif_list = results_dictionary[output_motif_csv]
    line_count = 0
    for line in output_motif_csv_open:
        line_clean = line.strip("\n")
        if line_clean in output_motif_list:
            line_count += 1
    count_list.append(line_count)
    all_motif_list = results_dictionary[output_all_motifs_csv]
    line_count = 0
    for line in output_all_motifs_csv_open:
        line_clean = line.strip("\n")
        if line_clean in all_motif_list:
            line_count += 1
    count_list.append(line_count)  
    # Checks to make sure all counts are in count_match.
    for count in count_list:
        position_count = count_list.index(count)
        assert count == count_match[position_count]        