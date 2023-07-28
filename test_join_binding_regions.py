"""test_join_binding_regions.py- Unit tests for
    join_binding_regions.py
"""
import pytest
# Imports all functions from join_binding_regions.py.
import join_binding_regions as j_b_r

"""
Tests add_events_to_dictionary() with different BED files.

Data Types:
    sample_1_bedfile -- BED file produced from Piranha or other CLIP
        read pileup tool. 
    region_dictionary -- Filled output dictionary with each read pileup
        region. Contains genomic location information as the value and
        the full bed line as the key. {"chrom_cord1_cord2": "line_clean"}.
    """

@pytest.mark.parametrize("sample_1_bedfile, region_dictionary", [
                         # Test 1: One line BED file
                         # sample_1_bedfile
                         (("test_files/join_binding_regions/test1.bed"),
                          # region_dictionary
                          ({"chr1_100_200_+": "chr1_100_200_X_15_+_0.00103521"})),
                         # Test 2: Multi-line BED file
                         # sample_1_bedfile
                         (("test_files/join_binding_regions/test2.bed"),
                          # region_dictionary
                          ({"chr1_100_200_+": "chr1_100_200_X_15_+_0.00103521",
                            "chr1_300_500_+": "chr1_300_500_X_15_+_0.00103521",
                            "chr1_600_1000_+": "chr1_600_1000_X_15_+_0.00103521",
                            "chr1_1500_2000_+": "chr1_1500_2000_X_15_+_0.00103521"}))
                         ])

def test_add_events_to_dictionary(sample_1_bedfile, region_dictionary):
    """
        GIVEN a BED file that needs to be added to a dictionary.
        WHEN the function adds each BED entry to a dictionary.
        THEN the dictionary is checked for the correct data in 
            the right format.
    """
    assert j_b_r.add_events_to_dictionary(sample_1_bedfile) == region_dictionary

"""
Tests output_file() with different BED files.

Data Types:
    region_dictionary_out -- Filled output dictionary with each read pileup
        region. Contains genomic location information as the value and
        the full bed line as the key. {"chrom_cord1_cord2": "line_clean"}  
    joined_gtf_file -- Output BED file. Contains nonrepetitive, comprehensive
        read pileup regions from each input sample
    results_list -- List of results expected in function produced joined_gtf_file. 
    count_match -- Expected count of matching lines in output and premade
        testing files.    
    """

@pytest.mark.parametrize("region_dictionary_out, joined_gtf_file, \
                         results_list, count_match", [
                         # Test 1: One line BED file.
                         # region_dictionary_out
                         (({"chr1_100_200_+": "chr1_100_200_X_15_+_0.00103521"}),
                          # joined_gtf_file
                          ("test_files/join_binding_regions/test1_out.bed"),
                          # results_list
                          (["chr1_pir_gene_100_200_15_+_0.00103521_gene_id \"chr1~100~200~+g\"; gene_version \"1\"; ",
                            "chr1_pir_transcript_100_200_15_+_0.00103521_gene_id \"chr1~100~200~+\"; gene_version \"1\"; transcript_id \"chr1~100~200~+t\"; "]),
                          #(["chr1_100_200_X_15_+_0.00103521"]),
                          # count_match  
                            (2)),
                         # Test 2: Multi-line BED file.
                         # region_dictionary_out
                         (({"chr1_100_200_+": "chr1_100_200_X_15_+_0.00103521",
                            "chr1_300_500_+": "chr1_300_500_X_15_+_0.00103521",
                            "chr1_600_1000_+": "chr1_600_1000_X_15_+_0.00103521",
                            "chr1_1500_2000_+": "chr1_1500_2000_X_15_+_0.00103521"}),
                          # joined_gtf_file
                          ("test_files/join_binding_regions/test2_out.bed"),
                          # results_list
                          (["chr1_pir_gene_100_200_15_+_0.00103521_gene_id \"chr1~100~200~+g\"; gene_version \"1\"; ",
                            "chr1_pir_transcript_100_200_15_+_0.00103521_gene_id \"chr1~100~200~+\"; gene_version \"1\"; transcript_id \"chr1~100~200~+t\"; ",
                            "chr1_pir_gene_300_500_15_+_0.00103521_gene_id \"chr1~300~500~+g\"; gene_version \"1\"; ",
                            "chr1_pir_transcript_300_500_15_+_0.00103521_gene_id \"chr1~300~500~+\"; gene_version \"1\"; transcript_id \"chr1~300~500~+t\"; ",
                            "chr1_pir_gene_600_1000_15_+_0.00103521_gene_id \"chr1~600~1000~+g\"; gene_version \"1\"; ",
                            "chr1_pir_transcript_600_1000_15_+_0.00103521_gene_id \"chr1~600~1000~+\"; gene_version \"1\"; transcript_id \"chr1~600~1000~+t\"; ",
                            "chr1_pir_gene_1500_2000_15_+_0.00103521_gene_id \"chr1~1500~2000~+g\"; gene_version \"1\"; ",
                            "chr1_pir_transcript_1500_2000_15_+_0.00103521_gene_id \"chr1~1500~2000~+\"; gene_version \"1\"; transcript_id \"chr1~1500~2000~+t\"; "]),
                          # count_match
                            (8))
                         ])

def test_output_file(region_dictionary_out, joined_gtf_file, 
                     results_list, count_match):
    """
        GIVEN a BED file dictionary.
        WHEN the function outputs the information in the dictionary
            to a BED output file
        THEN the output file is checked for the correct output by
            counting the number of perfectly matching lines.
    """
    # Sets line_count to 0. This is the count of how many
    # lines match in the output and the expected output file.
    line_count = 0
    # Executes function.
    j_b_r.output_file(region_dictionary_out, joined_gtf_file)
    # Opens testing output files to iterate through to check against
    # temparory output below.
    read1_opened = open(joined_gtf_file, 'r')
    # Checks each line of the output against the expected output
    # in the results_list. Counts should match.
    for line in read1_opened:
        # Makes line so it is easier to compare.
        line_clean = line.strip("\n").replace("\t","_")
        if line_clean in results_list:
            line_count += 1
    # Checks line counts.
    assert count_match == line_count 

"""
Tests join_binding_regions() with read pileup regions

Data Types:
    sample_1_bedfile -- BED file produced from Piranha or other CLIP
        read pileup tool.
    additional_bedfiles -- List of additional BED files (strings) 
        needed to join first BED file. BED files should be produced 
        from Piranha or other CLIP read pileup tool. Contains one 
        or more bed files.
    joined_gtf_file -- Output BED file. Contains nonrepetitive,
        comprehensive read pileup regions from each input sample.}.
    results_list -- List of results expected in function produced joined_gtf_file. 
    count_match -- Expected count of matching lines in output and premade
        testing files.    
    """

@pytest.mark.parametrize("sample_1_bedfile, additional_bedfiles, joined_gtf_file, \
                         results_list, count_match", [
                         # Test 1: Joined with one BED file, identical entries.
                         # sample_1_bedfile
                         (("test_files/join_binding_regions/test1.bed"),
                          # additional_bedfiles
                          (["test_files/join_binding_regions/test12.bed"]), 
                         # joined_gtf_file
                         ("test_files/join_binding_regions/test1_out.bed"),
                         # results_list
                         (["chr1_pir_gene_100_200_15_+_0.00103521_gene_id \"chr1~100~200~+g\"; gene_version \"1\"; ",
                            "chr1_pir_transcript_100_200_15_+_0.00103521_gene_id \"chr1~100~200~+\"; gene_version \"1\"; transcript_id \"chr1~100~200~+t\"; "]),
                         #(["chr1_100_200_X_15_+_0.00103521"]),
                         # count_match  
                         (2)),
                         # Test 2: Joined with two BED files, identical entries.
                         # sample_1_bedfile
                         (("test_files/join_binding_regions/test2.bed"),
                          # additional_bedfiles
                          (["test_files/join_binding_regions/test22.bed",
                            "test_files/join_binding_regions/test23.bed"]),
                          # joined_gtf_file
                          ("test_files/join_binding_regions/test2_out.bed"),
                          # results_list
                          (["chr1_pir_gene_100_200_15_+_0.00103521_gene_id \"chr1~100~200~+g\"; gene_version \"1\"; ",
                            "chr1_pir_transcript_100_200_15_+_0.00103521_gene_id \"chr1~100~200~+\"; gene_version \"1\"; transcript_id \"chr1~100~200~+t\"; ",
                            "chr1_pir_gene_300_500_15_+_0.00103521_gene_id \"chr1~300~500~+g\"; gene_version \"1\"; ",
                            "chr1_pir_transcript_300_500_15_+_0.00103521_gene_id \"chr1~300~500~+\"; gene_version \"1\"; transcript_id \"chr1~300~500~+t\"; ",
                            "chr1_pir_gene_600_1000_15_+_0.00103521_gene_id \"chr1~600~1000~+g\"; gene_version \"1\"; ",
                            "chr1_pir_transcript_600_1000_15_+_0.00103521_gene_id \"chr1~600~1000~+\"; gene_version \"1\"; transcript_id \"chr1~600~1000~+t\"; ",
                            "chr1_pir_gene_1500_2000_15_+_0.00103521_gene_id \"chr1~1500~2000~+g\"; gene_version \"1\"; ",
                            "chr1_pir_transcript_1500_2000_15_+_0.00103521_gene_id \"chr1~1500~2000~+\"; gene_version \"1\"; transcript_id \"chr1~1500~2000~+t\"; "]),
                          # count_match  
                          (8)),
                         # Test 3: Joined with two BED files, assorted overlapping entries with multiple
                         # chromosomes.
                         # sample_1_bedfile
                         (("test_files/join_binding_regions/test3.bed"),
                          # additional_bedfiles
                          (["test_files/join_binding_regions/test32.bed",
                            "test_files/join_binding_regions/test33.bed"]),
                          # joined_gtf_file
                          ("test_files/join_binding_regions/test3_out.bed"),
                          # results_list
                          (["chr1_pir_gene_100_250_15_+_0.00103521_gene_id \"chr1~100~250~+g\"; gene_version \"1\"; ",
                            "chr1_pir_transcript_100_250_15_+_0.00103521_gene_id \"chr1~100~250~+\"; gene_version \"1\"; transcript_id \"chr1~100~250~+t\"; ",
                            "chr1_pir_gene_100_1000_15_-_0.00111111_gene_id \"chr1~100~1000~-g\"; gene_version \"1\"; ",
                            "chr1_pir_transcript_100_1000_15_-_0.00111111_gene_id \"chr1~100~1000~-\"; gene_version \"1\"; transcript_id \"chr1~100~1000~-t\"; ",
                            "chr1_pir_gene_275_500_15_+_0.00103521_gene_id \"chr1~275~500~+g\"; gene_version \"1\"; ",
                            "chr1_pir_transcript_275_500_15_+_0.00103521_gene_id \"chr1~275~500~+\"; gene_version \"1\"; transcript_id \"chr1~275~500~+t\"; ",
                            "chr1_pir_gene_525_1000_15_+_0.00103521_gene_id \"chr1~525~1000~+g\"; gene_version \"1\"; ",
                            "chr1_pir_transcript_525_1000_15_+_0.00103521_gene_id \"chr1~525~1000~+\"; gene_version \"1\"; transcript_id \"chr1~525~1000~+t\"; ",
                            "chr1_pir_gene_1500_2200_15_+_0.00103521_gene_id \"chr1~1500~2200~+g\"; gene_version \"1\"; ",
                            "chr1_pir_transcript_1500_2200_15_+_0.00103521_gene_id \"chr1~1500~2200~+\"; gene_version \"1\"; transcript_id \"chr1~1500~2200~+t\"; ",                   
                            "chr1_pir_gene_2400_3100_15_+_0.00103521_gene_id \"chr1~2400~3100~+g\"; gene_version \"1\"; ",
                            "chr1_pir_transcript_2400_3100_15_+_0.00103521_gene_id \"chr1~2400~3100~+\"; gene_version \"1\"; transcript_id \"chr1~2400~3100~+t\"; ",
                            "chr1_pir_gene_3400_4100_15_+_0.00103521_gene_id \"chr1~3400~4100~+g\"; gene_version \"1\"; ",
                            "chr1_pir_transcript_3400_4100_15_+_0.00103521_gene_id \"chr1~3400~4100~+\"; gene_version \"1\"; transcript_id \"chr1~3400~4100~+t\"; ",
                            "chr1_pir_gene_4400_5100_15_+_0.00103521_gene_id \"chr1~4400~5100~+g\"; gene_version \"1\"; ",
                            "chr1_pir_transcript_4400_5100_15_+_0.00103521_gene_id \"chr1~4400~5100~+\"; gene_version \"1\"; transcript_id \"chr1~4400~5100~+t\"; ",
                            "chr1_pir_gene_5500_6000_15_+_0.00103521_gene_id \"chr1~5500~6000~+g\"; gene_version \"1\"; ",
                            "chr1_pir_transcript_5500_6000_15_+_0.00103521_gene_id \"chr1~5500~6000~+\"; gene_version \"1\"; transcript_id \"chr1~5500~6000~+t\"; ",
                            "chr2_pir_gene_50_650_15_+_0.00111111_gene_id \"chr2~50~650~+g\"; gene_version \"1\"; ",
                            "chr2_pir_transcript_50_650_15_+_0.00111111_gene_id \"chr2~50~650~+\"; gene_version \"1\"; transcript_id \"chr2~50~650~+t\"; ",
                            "chr5_pir_gene_100_1000_15_+_0.00111111_gene_id \"chr5~100~1000~+g\"; gene_version \"1\"; ",
                            "chr5_pir_transcript_100_1000_15_+_0.00111111_gene_id \"chr5~100~1000~+\"; gene_version \"1\"; transcript_id \"chr5~100~1000~+t\"; ",
                            "chr6_pir_gene_100_700_15_+_0.00111111_gene_id \"chr6~100~700~+g\"; gene_version \"1\"; ",
                            "chr6_pir_transcript_100_700_15_+_0.00111111_gene_id \"chr6~100~700~+\"; gene_version \"1\"; transcript_id \"chr6~100~700~+t\"; ",
                            "chr6_pir_gene_800_1000_15_+_0.00111111_gene_id \"chr6~800~1000~+g\"; gene_version \"1\"; ",
                            "chr6_pir_transcript_800_1000_15_+_0.00111111_gene_id \"chr6~800~1000~+\"; gene_version \"1\"; transcript_id \"chr6~800~1000~+t\"; ",                          
                            "chr6_pir_gene_100_6000_15_-_0.00111111_gene_id \"chr6~100~6000~-g\"; gene_version \"1\"; ",
                            "chr6_pir_transcript_100_6000_15_-_0.00111111_gene_id \"chr6~100~6000~-\"; gene_version \"1\"; transcript_id \"chr6~100~6000~-t\"; "]),
                          # count_match  
                          (28)),
                         ])

def test_join_binding_regions(sample_1_bedfile, additional_bedfiles, joined_gtf_file,
                              results_list, count_match):
    """
        GIVEN a BED file that needs to be joined with additional BED files.
        WHEN the function checks each annotation for overlapping entries,
            removes overlapping entries, and then saves the data in a dictionary.
        THEN the output file is checked for the correct output by
            counting the number of perfectly matching lines.
    """
    # Sets line_count to 0. This is the count of how many
    # lines match in the output and the expected output file.
    line_count = 0
    # Executes function.
    j_b_r.join_binding_regions(sample_1_bedfile, additional_bedfiles, joined_gtf_file)
    # Opens testing output files to iterate through to check against
    # temporary output below.
    read1_opened = open(joined_gtf_file, 'r')
    # Checks each line of the output against the expected output
    # in the results_list. Counts should match.
    for line in read1_opened:
        # Makes line so it is easier to compare.
        line_clean = line.strip("\n").replace("\t","_")
        if line_clean in results_list:
            line_count += 1
    # Checks line counts.
    assert count_match == line_count
