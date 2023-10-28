"""test_make_random_bed_for_motif.py- Unit tests for
    make_random_bed_for_motif.py
"""
import pytest
# Imports all functions from make_random_bed_for_motif.py.
import make_random_bed_for_motif as m_r_b_m

"""
Tests get_exons_from_gtf using assorted GTF files.

Data Types:
    gtf_file -- Gencode/Ensemlb or other GTF file
        for the organism of interest. Should annotate
        the FASTA file used for alignment.
    exon_list -- List of exons from GTF file.
    results_list -- List of expected results from function.
    count_match -- List of expected count of matching lines in output 
        and premade testing files.    
"""

@pytest.mark.parametrize("gtf_file, results_list, count_match", [
                         # Test 1: Single gene GTF with one exon.
                         # gtf_file
                         (("test_files/make_random_bed_for_motif/gtf_test1.tsv"),
                          # results_list
                          (["chr1_251_350_+"]),
                          # count_match
                          (1)),
                         # Test 2: Multiple gene GTF with multiple exons.
                         # gtf_file
                         (("test_files/make_random_bed_for_motif/gtf_test2.tsv"),
                          # results_list
                          (["chr1_251_350_+", "chr3_300_500_+"]),
                          # count_match  
                          (2)),
                         # Test 3: Multiple gene GTF with multiple exons, including
                         # negative strand, with one exon below 250.
                         # gtf_file
                         (("test_files/make_random_bed_for_motif/gtf_test3.tsv"),
                          # results_list
                          (["chr1_255_260_-", "chr2_4500_4600_+"]),
                          # count_match  
                          (2)),
                         ])

def test_get_exons_from_gtf(gtf_file, results_list, count_match):
    """
        GIVEN a GTF file.
        WHEN the function pulls out exons from the GTF file and outputs
            them in a list.
        THEN the output list is checked for the correct entries in
            the correct number.
    """
    # Count used to compare to count_match.
    line_count = 0
    # Executes function.
    exon_list = m_r_b_m.get_exons_from_gtf(gtf_file)
    # Iterates through list and checks number of correct
    # entries versus results_list.
    for exon in exon_list:
        if (exon in results_list):
            line_count += 1
    assert count_match == line_count

"""
Tests generate_regions using assorted GTF files.

Data Types:
    gtf_file -- Gencode/Ensemlb or other GTF file
        for the organism of interest. Should annotate
        the FASTA file used for alignment.
    output_file_prefix -- Prefix for output files.
    random_regions -- Number of random regions to generate.
    replicates -- Number of BED files to generate.
    results_list -- List of expected results from function.
    count_match -- List of expected count of matching lines in output 
        and premade testing files.    
"""

@pytest.mark.parametrize("gtf_file, output_file_prefix, random_regions, replicates, \
                          results_list, count_match_list", [
                         # Test 1: Single gene GTF with one exon.
                         # gtf_file
                         (("test_files/make_random_bed_for_motif/gtf_test_rand1.tsv"),
                          # output_file_prefix
                          ("test_files/make_random_bed_for_motif/test1"),
                          # random_regions
                          (1),
                          # replicates
                          (1),
                         # results_list
                         (["chr1_251_350_+"]),
                         # count_match_list 
                         ([1])),
                         # Test 2: Multiple gene GTF with multiple exons.
                         # gtf_file
                         (("test_files/make_random_bed_for_motif/gtf_test_rand2.tsv"),
                          # output_file_prefix
                          ("test_files/make_random_bed_for_motif/test2"),
                          # random_regions
                          (3),
                          # replicates
                          (1),
                          # results_list
                          (["chr1_251_350_+", "chr3_300_500_+"]),
                          # count_match_list 
                          ([3])),
                          # Test 3: Multiple gene GTF with multiple exons, including
                          # negative strand. One exon is small (< 250) and should be 
                          # filtered out.
                          # gtf_file
                          (("test_files/make_random_bed_for_motif/gtf_test_rand3.tsv"),
                            # output_file_prefix
                            ("test_files/make_random_bed_for_motif/test3"),
                            # random_regions
                            (4),
                            # replicates
                            (1),
                            # results_list
                            (["chr1_500_501_+", "chr1_1000_1500_-", "chr1_2000_2001_+", "chr2_4500_4600_+"]),
                            # count_match_list  
                            ([4])),
                          # Test 4: Single gene GTF with one exon and two replicates.
                          # gtf_file
                          (("test_files/make_random_bed_for_motif/gtf_test_rand1.tsv"),
                           # output_file_prefix
                           ("test_files/make_random_bed_for_motif/test4"),
                           # random_regions
                           (1),
                           # replicates
                           (2),
                           # results_list
                           (["chr1_251_350_+"]),
                           # count_match_list 
                           ([1, 1])),
                          # Test 5: Multiple gene GTF with multiple exons, including
                          # negative strand. One exon is small (< 250) and should be 
                          # filtered out.
                          # gtf_file
                          (("test_files/make_random_bed_for_motif/gtf_test_rand3.tsv"),
                           # output_file_prefix
                           ("test_files/make_random_bed_for_motif/test5"),
                           # random_regions
                           (4),
                           # replicates
                           (4),
                           # results_list
                           (["chr1_500_501_+", "chr1_1000_1500_-", "chr1_2000_2001_+", "chr2_4500_4600_+"]),
                           # count_match_list 
                           ([4, 4, 4, 4])),
                        ])

def test_generate_regions(gtf_file, output_file_prefix, random_regions, replicates,
                          results_list, count_match_list):
    """
        GIVEN a GTF file and a need to generate random regions.
        WHEN the function generates random regions and outputs
            them in a BED file.
        THEN the output BED file is checked for the correct entries in
            the correct number for each desired replicate.
    """
    # Used to count number of replicates tested.
    replicate_count = 0
    count = 0
    # List used to keep track of the correct line_count matches for 
    # each replicate.
    line_count_list = []
    # Executes function.
    m_r_b_m.generate_regions(gtf_file, random_regions, replicates, output_file_prefix)
    # Iterates through list and looks for exon entries based on
    # the random region algorythm parameters. Chromosome and strand
    # must match and one of the coordinates must be with 250 nucleotides.
    # For testing, exons should be spaced out so that they are not counted
    # twice or on a different chromosome or strand.
    while (replicate_count < replicates):
        bed_file = open(output_file_prefix + "_random_" + str(replicate_count) + ".bed", 'r')
        # Count used to compare to count_match.
        line_count = 0
        for line in bed_file:
            line_list = line.strip("\n").split("\t")
            for result in results_list:
                result_split = result.split("_")
                if ((result_split[0] == line_list[0]) and (result_split[3] == line_list[5])):
                    if ((abs(int(result_split[1]) - int(line_list[1])) < 251)
                        or (abs(int(result_split[2]) - int(line_list[2])) < 251)):
                        line_count += 1
        line_count_list.append(line_count)        
        replicate_count += 1        
    # Count used to grab correct items in list for assertions.
    list_count = 0
    # Checks to make sure all counts are in count_match_list. Use different
    # number of counts for each replicate to ensure all counts match.            
    for count in line_count_list:
        assert count == count_match_list[list_count]
        list_count += 1   

                              
   

