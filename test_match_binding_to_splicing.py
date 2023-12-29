"""test_match_binding_to_splicing.py- Unit tests for
    match_binding_to_splicing.py
"""
import pytest
# Imports all functions from match_binding_to_splicing.py.
import match_binding_to_splicing as m_b_t_s

"""
Tests collect_regions() with different numbers of binding regions.

Data Types:
    region_csv -- CSV file of random motif enrichment values
        for each binding region from all input files. Was generated
        from the random_enrichment_motify.py script. File name-
        output_file_prefix + "_motifs_per_region.csv". Format -
        Region Name, Random Motif Enrichment Value, Normalized. 
    region_dictionary -- Dictionary of binding regions with
        region name as key and motif values as value.
        {region_name: [motif enrichment, normalized]}     
    """

@pytest.mark.parametrize("region_csv, region_dictionary", [
                         # Test 1: One binding region.
                         # region_csv
                         (('test_files/match_binding_to_splicing/test1_motifs_per_region.csv'),
                           # region_dictionary
                           ({'name1': ['1.0', 'Yes']})),
                         # Test 2: Multiple binding regions.
                         # region_csv
                         (('test_files/match_binding_to_splicing/test2_motifs_per_region.csv'),
                           # region_dictionary
                           ({'name1': ['0.25', 'No'],
                             'name2': ['0.25', 'No'],
                             'name3': ['0.25', 'No'],
                             'name41': ['0.25', 'No'],
                             'name42': ['0.5', 'No'],
                             'name43': ['0.5', 'No']})),          
                         ])

def test_collect_regions(region_csv, region_dictionary):
    """
        GIVEN a random enrichment motif value file.
        WHEN the function parses the file and adds the motif value information
            to a dictionary.
        THEN the output dictionary is checked for the correct motif values.
    """
    assert m_b_t_s.collect_regions(region_csv) == (region_dictionary)

"""
Tests make_splicing_dictionary() with different numbers of splicing events.

Data Types:
    rmats _file -- rMATS alternative skipped exon splicing 
        output file (SE.MATS.JCEC.txt file). Should have 
        been generated using the same target organism 
        GTF file. All other types of alternative splicing 
        can be input except for retained intron (different 
        number of columns).
    output -- String prefix for output file. Should contain
        the full path to the output directory if a different
        directory is being used.    
    splicing_dictionary -- Dictionary of splicing events
        with exon coordinates as keys and additional
        splicing information as the value.
        {"gene-chromosome-exon_start-exon_end-strand":
         [ensemble_id, up-start, up-end, down-start, down-end, 
          delta_psi, regulation]}
    filtered_splicing_out -- rMATS style output file that is
        filtered by delta PSI and an FDR of 10% or less.
        File name- output + "filtered_SE.MATS.JCEC.txt".
    results_check -- Text file with results expected in function 
        produced files.    
    count_match -- Expected count of matching lines in output 
        and premade testing files.        
    """

@pytest.mark.parametrize("rmats_file, output, splicing_dictionary, filtered_splicing_out, \
                          results_check, count_match", [
                         # Test 1: One splicing event that does not pass.
                         # rmats_file
                         (('test_files/match_binding_to_splicing/test1_SE.MATS.JCEC.txt'),
                          # output
                          ('test_files/match_binding_to_splicing/test1_'),
                          # splicing_dictionary
                          ({}),
                          # filtered_splicing_out
                          ('test_files/match_binding_to_splicing/test1_filtered_SE.MATS.JCEC.txt'),
                          # results_check
                          ('test_files/match_binding_to_splicing/test1_SE.MATS.JCEC.check.txt'),
                          # count_match
                          (1)),
                          # Test 2: Three splicing events, one that has no passing characteristics,
                          # on that is significant but not high enough delta PSI, and one that
                          # is significant and high enough delta PSI.
                          # rmats_file
                          (('test_files/match_binding_to_splicing/test2_SE.MATS.JCEC.txt'),
                           # output
                           ('test_files/match_binding_to_splicing/test2_'),
                           # splicing_dictionary
                           ({'Efemp2~chr19~5525179~5525206~+': ['ENSMUSG00000024909.16',
                                                                '5524774', '5524802',
                                                                '5525435', '5525484',
                                                                '0.11', 'Exclusion']}),
                           # filtered_splicing_out
                           ('test_files/match_binding_to_splicing/test2_filtered_SE.MATS.JCEC.txt'),
                           # results_check
                           ('test_files/match_binding_to_splicing/test2_SE.MATS.JCEC.check.txt'),
                           # count_match
                           (2)),
                          # Test 3: Five splicing events, one that has no passing characteristics,
                          # on that is significant but not high enough delta PSI, and three that
                          # are significant, have a high enough delta PSI, and a mix of inclusion
                          # and exclusion.
                          # rmats_file
                          (('test_files/match_binding_to_splicing/test3_SE.MATS.JCEC.txt'),
                           # output
                           ('test_files/match_binding_to_splicing/test3_'),
                           # splicing_dictionary
                           ({'Efemp2~chr19~5525179~5525206~+': ['ENSMUSG00000024909.16',
                                                                '5524774', '5524802',
                                                                '5525435', '5525484',
                                                                '0.11', 'Exclusion'],
                             'Efemp21~chr19~5525088~5525206~+': ['ENSMUSG00000024909.17',
                                                                 '5524715', '5524802',
                                                                 '5525435', '5525484',
                                                                 '-0.51', 'Inclusion'],
                             'Efemp22~chr19~5525179~5525206~+': ['ENSMUSG00000024909.18',
                                                                 '5524774', '5524802',
                                                                 '5525435', '5525484',
                                                                 '0.60', 'Exclusion']}),
                           # filtered_splicing_out
                           ('test_files/match_binding_to_splicing/test3_filtered_SE.MATS.JCEC.txt'),
                           # results_check
                           ('test_files/match_binding_to_splicing/test3_SE.MATS.JCEC.check.txt'),
                           # count_match
                           (4)),         
                         ])

def test_make_splicing_dictionary(rmats_file, output, splicing_dictionary, filtered_splicing_out,
                                  results_check, count_match):
    """
        GIVEN a random enrichment motif value file.
        WHEN the function parses the file and adds the motif value information
            to a dictionary.
        THEN the output dictionary is checked for the correct motif values.
    """
    # Counts number of lines in the output file.
    line_count = 0
    assert m_b_t_s.make_splicing_dictionary(rmats_file, output) == (splicing_dictionary)
    out1_opened = open(filtered_splicing_out, 'r')
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
    # Checks line counts.
    assert count_match == line_count
    out1_opened.close()
    check1_opened.close()

"""
Tests match_binding_to_splicing() with different numbers of 
    splicing and binding events.

Data Types:
    binding_regions -- DESeq2 output file with binding
        regions. Should be annotated and filtered file,
        although any DESeq2 output file can be used.
    splicing_dictionary -- Dictionary of splicing events
        with exon coordinates as keys and additional
        splicing information as the value.
        {"gene-chromosome-exon_start-exon_end-strand":
         [ensemble_id, up-start, up-end, down-start, down-end, 
          delta_psi, regulation]}
    distance -- Integer distance skipped exon outer
        cordinates should be to be considered within
        regulated distance of binding region.
    output -- String prefix for output file. Should contain
        the full path to the output directory if a different
        directory is being used.
    region_csv -- CSV file of random motif enrichment values
        for each binding region from all input files. Was generated
        from the random_enrichment_motify.py script. File name-
        output_file_prefix + "_motifs_per_region.csv". Format -
        Region Name, Random Motif Enrichment Value, Normalized.     
    splicing_binding_out -- CSV file with splicing and binding
        information for each region.    
    results_check -- Text file with results expected in function 
        produced files.    
    count_match -- Expected count of matching lines in output 
        and premade testing files.        
    """

@pytest.mark.parametrize("binding_regions, splicing_dictionary, distance, \
                          output, region_csv, splicing_binding_out, \
                          results_check, count_match", [
                         # Test 1: No splicing events.
                         # binding_regions
                         (('test_files/match_binding_to_splicing/test1_annotated_binding_regions.csv'),
                          # splicing_dictionary
                          ({}),
                          # distance
                          (100),
                          # output
                          ('test_files/match_binding_to_splicing/test1_'),
                          # region_csv
                          ('test_files/match_binding_to_splicing/test1_motifs_per_region.csv'),
                          # splicing_binding_out
                          ('test_files/match_binding_to_splicing/test1_splicing_binding.csv'),
                          # results_check
                          ("test_files/match_binding_to_splicing/test1_splicing_binding_check.csv"),
                          # count_match
                          (1)),
                         # Test 2: One matching splicing event with a motif value.
                         # binding_regions
                         (('test_files/match_binding_to_splicing/test1_annotated_binding_regions.csv'),
                          # splicing_dictionary
                          ({'Efemp2~chr1~399~451~+': ['ENSMUSG00000024909.16',
                                                      '5524774', '5524802',
                                                      '5525435', '5525484',
                                                      '0.11', 'Exclusion']}),
                          # distance
                          (100),
                          # output
                          ('test_files/match_binding_to_splicing/test2_'),
                          # region_csv
                          ('test_files/match_binding_to_splicing/test1_motifs_per_region.csv'),
                          # splicing_binding_out
                          ('test_files/match_binding_to_splicing/test2_splicing_binding.csv'),
                          # results_check
                          ("test_files/match_binding_to_splicing/test2_splicing_binding_check.csv"),
                          # count_match
                          (2)),
                         # Test 3: Multiple matching splicing event with and without motif values. The
                         # events also test multiple configurations of canonical and non-canonical
                         # splicing events. See below for more details.
                         # binding_regions
                         (('test_files/match_binding_to_splicing/test3_annotated_binding_regions.csv'),
                          # splicing_dictionary
                          # Tests upstream and canonical.
                          ({'Efemp2~chr1~399~451~+': ['ENSMUSG00000024909.16',
                                                      '5524774', '5524802',
                                                      '5525435', '5525484',
                                                      '0.11', 'Exclusion'],
                             # Tests upstream and non-canonical.                         
                             'Efemp2~chr2~399~451~-': ['ENSMUSG00000024909.17',
                                                        '5524715', '5524802',
                                                        '5525435', '5525484',
                                                        '-0.51', 'Inclusion'],
                             # Tests inside (from downstream) and canonical.                           
                             'Efemp2~chr3~150~301~+': ['ENSMUSG00000024909.18',
                                                        '5524774', '5524802',
                                                        '5525435', '5525484',
                                                        '0.60', 'Exclusion'],
                             # Tests inside (from upstream) and non-canonical.                       
                             'Efemp2~chr4~75~125~-': ['ENSMUSG00000024909.18',
                                                      '5524774', '5524802',
                                                      '5525435', '5525484',
                                                      '-0.60', 'Inclusion'],
                             # Tests downstream and canonical.                         
                             'Efemp2~chr5~50~75~+': ['ENSMUSG00000024909.18',
                                                       '5524774', '5524802',
                                                       '5525435', '5525484',
                                                       '-0.60', 'Inclusion'],
                             # Tests downstream and non-canonical.                                                 
                             'Efemp2~chr5~25~50~+': ['ENSMUSG00000024909.18',
                                                       '5524774', '5524802',
                                                       '5525435', '5525484',
                                                       '0.60', 'Exclusion'],
                             # Tests distance.                          
                             'Efemp2~chr5~401~501~+': ['ENSMUSG00000024909.18',
                                                       '5524774', '5524802',
                                                       '5525435', '5525484',
                                                       '0.60', 'Exclusion'],
                             # Tests matching inside coordinates.                          
                             'Efemp2~chr7~100~300~+': ['ENSMUSG00000024909.18',
                                                       '5524774', '5524802',
                                                       '5525435', '5525484',
                                                       '0.60', 'Exclusion'],
                             # Tests matching upstream coordinate.                          
                             'Efemp2~chr7~50~100~+': ['ENSMUSG00000024909.18',
                                                       '5524774', '5524802',
                                                       '5525435', '5525484',
                                                       '0.60', 'Exclusion'],
                             # Tests matching downstream coordinate and non-canonical.                          
                             'Efemp2~chr7~300~350~+': ['ENSMUSG00000024909.18',
                                                       '5524774', '5524802',
                                                       '5525435', '5525484',
                                                       '-0.60', 'Inclusion'],
                             # Tests non-matching strand.                                                                                                                                         
                             'Efemp2~chr5~75~125~-': ['ENSMUSG00000024909.18',
                                                      '5524774', '5524802',
                                                      '5525435', '5525484',
                                                      '-0.60', 'Inclusion'],
                             # Tests non-matching chromosome.                         
                             'Efemp2~chr6~75~125~-': ['ENSMUSG00000024909.18',
                                                      '5524774', '5524802',
                                                      '5525435', '5525484',
                                                      '-0.60', 'Inclusion']}),
                          # distance
                          (100),
                          # output
                          ('test_files/match_binding_to_splicing/test3_'),
                          # region_csv
                          ('test_files/match_binding_to_splicing/test3_motifs_per_region.csv'),
                          # splicing_binding_out
                          ('test_files/match_binding_to_splicing/test3_splicing_binding.csv'),
                          # results_check
                          ("test_files/match_binding_to_splicing/test3_splicing_binding_check.csv"),
                          # count_match
                          (10)),           
                         ])

def test_match_binding_to_splicing(binding_regions, splicing_dictionary, distance,
                                  output, region_csv, splicing_binding_out, 
                                  results_check, count_match):
    """
        GIVEN a random enrichment motif value file.
        WHEN the function parses the file and adds the motif value information
            to a dictionary.
        THEN the output file is checked for the correct motif values.
    """
    # Counts number of lines in the output file.
    line_count = 0
    m_b_t_s.match_binding_to_splicing(binding_regions, splicing_dictionary,
                                      distance, output, region_csv)
    out1_opened = open(splicing_binding_out, 'r')
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
    # Checks line counts.
    assert count_match == line_count
    out1_opened.close()
    check1_opened.close()    