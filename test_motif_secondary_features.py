"""test_motif_secondary_features.py- Unit tests for
    motif_secondary_features.py
"""
import pytest
# Imports all functions from make_random_bed_for_motif.py.
import motif_secondary_features as m_s_f

"""
Tests organize_nucleotides with various nucleotide lists.

    Data Types:
    fasta_file_name -- String, name of FASTA file being analyzed
    kmer_size -- Size of kmer to be analyzed.
    motif_list -- List of motifs to be analyzed. Default is YGCY.
    motif_title -- Title for motif output file. Default is YGCY.
    region_dictionary -- Dictionary of nucleotide lists for the nucleotides
        surrounding each motif and the motifs themselves. The data is arranged
        in a list of four lists. The descriptions of each list are listed below.
        "nucleotide_list"- List, string list of nucleotides, seven nucleotides
            around all motifs.
        "nucleotide_intra_list"- List, string list of nucleotides between motifs
            in a motif group. A motif group is defined as motifs within seven 
            nucleotides of each other.
        "nucleotide_inter_list"- List, string list of nucleotides between motifs
            separated by more than seven nucleotides.
        "motifs_group_list"- List, string list of motifs and surrounding 
            nucleotides within a group.
        region_size -- Integer, size of binding region (last variable).    
        {"region_name": [nucleotide_list, nucleotide_intra_list,
                         nucleotide_inter_list, motifs_group_list,
                         region_size].     
"""

@pytest.mark.parametrize("fasta_file_name, kmer_size, motif_list, \
                          motif_title, region_dictionary", [
                        # Test 1: Single region with two motif groups. One
                        # group with two motifs, the other with one. Also
                        # tests a motif coming up first and last in the
                        # sequence.
                        # fasta_file_name
                        (("test_files/motif_secondary_features/test1.fasta"),
                        # kmer_size
                        (4),
                        # motif_list
                        (["CGCC", "CGCT", "TGCC", "TGCT"]),
                        # motif_title
                        ("YGCY"),
                        # region_dictionary
                        ({"name1": [['AAAATTTTTTTTTTTTTT'], 
                                    ['AAAA'], ['TTTTTTTT'], 
                                    ['TGCTAAAACGCCTTTTTTT', 'TTTTTTTCGCC'], 
                                    24]})),
                        # Test 2: Single region with two motif groups. One
                        # group with two motifs, the other with a double
                        # motif. Also checks for motifs that are not at
                        # beginning or end but the surrounding nucleotides
                        # would be.
                        # fasta_file_name
                        (("test_files/motif_secondary_features/test2.fasta"),
                        # kmer_size
                        (4),
                        # motif_list
                        (["CGCC", "CGCT", "TGCC", "TGCT"]),
                        # motif_title
                        ("YGCY"),
                        # region_dictionary
                        ({"name1": [['CCCAAAATTTTTTTTTTTTTTGGGG'], 
                                    ['AAAA'], ['TTTTTTTT'], 
                                    ['CCCTGCTAAAACGCCTTTTTTT', 
                                     'TTTTTTTCGCCCGCCGGGG'], 
                                    35]})),
                        # Test 3: Single region with three motif groups. One
                        # group with two motifs, one group with four, 
                        # the other with a double motif. Also checks for 
                        # motifs that are not at beginning or end but the 
                        # surrounding nucleotides would be.
                        # fasta_file_name
                        (("test_files/motif_secondary_features/test3.fasta"),
                        # kmer_size
                        (4),
                        # motif_list
                        (["CGCC", "CGCT", "TGCC", "TGCT"]),
                        # motif_title
                        ("YGCY"),
                        # region_dictionary
                        ({"name1": [['CCCCTTTTTTAGTGTGTG'], 
                                    ['TTTTTTA'], [], 
                                    ['CCCCCGCCTTTTTTCGCTATGCCGTGTGTG'], 
                                    33]})),
                        # Test 4: Single region with two motif groups. One
                        # group with two double motifs, the other with a 
                        # double motif and a single space then another motif. 
                        # All other surrounding regions are larger than seven-
                        # beginning, middle, and the end.
                        # fasta_file_name
                        (("test_files/motif_secondary_features/test4.fasta"),
                        # kmer_size
                        (4),
                        # motif_list
                        (["CGCC", "CGCT", "TGCC", "TGCT"]),
                        # motif_title
                        ("YGCY"),
                        # region_dictionary
                        ({"name1": [['CCCCCCCAAAAAAAAAAAATTAGGGGGGG'], 
                                    ['A'], ['AAAAAAAAATT'], 
                                    ['CCCCCCCTGCTGCCAAAAAAA',
                                     'AAAAATTCGCTGCTATGCTGGGGGGG'], 
                                    52]})),
                        # Test 5: Single region with two motif groups. One
                        # group with a double motif with three motifs, and
                        # multiple other spaced out motifs. Then a single
                        # motif in a closely followed group.
                        # fasta_file_name
                        (("test_files/motif_secondary_features/test5.fasta"),
                        # kmer_size
                        (4),
                        # motif_list
                        (["CGCC", "CGCT", "TGCC", "TGCT"]),
                        # motif_title
                        ("YGCY"),
                        # region_dictionary
                        ({"name1": [['AAGGGACCCCCCCCCCCCCC'], 
                                    ['AAGGGA'], ['CCCCCCCC'], 
                                    ['CGCCGCTGCTAATGCTGGGTGCTATGCTCCCCCCC',
                                     'CCCCCCCTGCT'], 
                                    40]})),
                        # Test 6: Multiple binding regions with assorted
                        # nucleotide patterns. Taken from test 1,3,5. One
                        # region also has no motifs.
                        # fasta_file_name
                        (("test_files/motif_secondary_features/test6.fasta"),
                        # kmer_size
                        (4),
                        # motif_list
                        (["CGCC", "CGCT", "TGCC", "TGCT"]),
                        # motif_title
                        ("YGCY"),
                        # region_dictionary
                        ({"name1": [['AAAATTTTTTTTTTTTTT'], 
                                    ['AAAA'], ['TTTTTTTT'], 
                                    ['TGCTAAAACGCCTTTTTTT', 'TTTTTTTCGCC'], 
                                    24],
                          "name3": [['CCCCTTTTTTAGTGTGTG'], 
                                    ['TTTTTTA'], [], 
                                    ['CCCCCGCCTTTTTTCGCTATGCCGTGTGTG'], 
                                    33],         
                          "name5": [['AAGGGACCCCCCCCCCCCCC'], 
                                    ['AAGGGA'], ['CCCCCCCC'], 
                                    ['CGCCGCTGCTAATGCTGGGTGCTATGCTCCCCCCC',
                                     'CCCCCCCTGCT'], 
                                    40],
                          "name6": [[], [], [], [], 9]})),          
                        # Test 7: Multiple lines with no motifs.
                        # fasta_file_name
                        (("test_files/motif_secondary_features/test7.fasta"),
                        # kmer_size
                        (4),
                        # motif_list
                        (["CGCC", "CGCT", "TGCC", "TGCT"]),
                        # motif_title
                        ("YGCY"),
                        # region_dictionary
                        ({"name1": [[], [], [], [], 22],
                          "name6": [[], [], [], [], 9]})),                                                                         
                        ])

def test_organize_nucleotides(fasta_file_name, kmer_size, motif_list,
                              motif_title, region_dictionary):
    """
        GIVEN list of motifs, a given motif size, a motif title, and a
            FASTA file of RBP binding regions that need to be analyzed.
        WHEN the function organizes the nucleotides surrounding each motif
            and the motifs themselves into four lists.
        THEN the output lists are checked for the correct nucleotides and
            motifs.
    """
    # Executes function and checks output.
    assert m_s_f.organize_nucleotides(fasta_file_name, kmer_size, motif_list,
                                      motif_title) == (region_dictionary, 
                                                       motif_title, motif_list)

"""
Tests get_motif_counts with various nucleotide lists.

    Data Types:
    region_dictionary -- Dictionary of nucleotide lists for the nucleotides
        surrounding each motif and the motifs themselves. The data is arranged
        is a list of four lists and an integer. The descriptions of each data
        structure is listed below. See organize_nucleotides for details.
        "nucleotide_list"- List, string list of nucleotides, seven nucleotides
            around all motifs.
        "nucleotide_intra_list"- List, string list of nucleotides between motifs
            in a motif group. A motif group is defined as motifs within seven 
            nucleotides of each other.
        "nucleotide_inter_list"- List, string list of nucleotides between motifs
            separated by more than seven nucleotides.
        "motifs_group_list"- List, string list of motifs and surrounding 
            nucleotides within a group.
        "region_size"- Integer, size of binding region (last variable).    
        {"region_name": [nucleotide_list, nucleotide_intra_list,
                         nucleotide_inter_list, motifs_group_list,
                         region_size].
        motif_list -- List of motifs to be analyzed. Default is YGCY.
        region_stats_dictionary -- Dictionary with motif and secondary
            statistics for each region. Dictionary is organized by region.
            See nucleotide_content and motif_content for details.
            {"region_name": [nucleotide_list, nucleotide_intra_list,
                             nucleotide_inter_list, motifs_group_list,
                             region_size]}. 
        file_stats_dictionary -- Dictionary with motif and seconardary
            statistics for each input file. See aggregate_file_stats
            for details.
            {"file_name": [nucleotide_count, nucleotide_intra_count,
                           nucleotide_inter_count, motif_region_counts]}.  
"""

@pytest.mark.parametrize("region_dictionary, motif_list, fasta_file_name, \
                         region_stats_dictionary, file_stats_dictionary", [
                        # Test 1: Single region with two motif groups. One
                        # group with two motifs, the other with one. Also
                        # tests a motif coming up first and last in the
                        # sequence.
                        # region_dictionary
                        (({"name1": [['AAAATTTTTTTTTTTTTT'], 
                                    ['AAAA'], ['TTTTTTTT'], 
                                    ['TGCTAAAACGCCTTTTTTT', 'TTTTTTTCGCC'], 
                                    24]}), 
                        # motif_list
                        (["CGCC", "CGCT", "TGCC", "TGCT"]),
                        # fasta_file_name
                        ("test_files/motif_secondary_features/test1.fasta"),
                        # region_stats_dictionary
                        ({'name1': [{'A': 4, 'C': 0, 'G': 0, 'T': 14, 'Total': 18, 
                                     'A_percentage': 22.22222222222222, 'C_percentage': 0.0, 
                                     'G_percentage': 0.0, 'T_percentage': 77.77777777777779, 
                                     'Average_length': 18.0}, 
                                     {'A': 4, 'C': 0, 'G': 0, 'T': 0, 'Total': 4, 
                                      'A_percentage': 100.0, 'C_percentage': 0.0, 
                                      'G_percentage': 0.0, 'T_percentage': 0.0, 
                                      'Average_length': 4.0}, 
                                      {'A': 0, 'C': 0, 'G': 0, 'T': 8, 'Total': 8, 
                                       'A_percentage': 0.0, 'C_percentage': 0.0, 
                                       'G_percentage': 0.0, 'T_percentage': 100.0, 
                                       'Average_length': 8.0}, 
                                       {'CGCC': 2, 'CGCT': 0, 'TGCC': 0, 'TGCT': 1, 
                                        'Total_motif': 3, 'A': 4, 'C': 0, 'G': 0, 'T': 14, 
                                        'Total_nucleotide': 18, 
                                        'A_rate_per_motif_per_group': 0.14814814814814814, 
                                        'A_per_group': 2.0, 'C_rate_per_motif_per_group': 0.0, 
                                        'C_per_group': 0.0, 'G_rate_per_motif_per_group': 0.0, 
                                        'G_per_group': 0.0, 'T_rate_per_motif_per_group': 0.5185185185185185, 
                                        'T_per_group': 7.0, 'A_per_CGCC': 2.0, 'C_per_CGCC': 0.0, 
                                        'G_per_CGCC': 0.0, 'T_per_CGCC': 7.0, 'A_per_CGCT': 0, 
                                        'C_per_CGCT': 0, 'G_per_CGCT': 0, 'T_per_CGCT': 0, 
                                        'A_per_TGCC': 0, 'C_per_TGCC': 0, 'G_per_TGCC': 0, 
                                        'T_per_TGCC': 0, 'A_per_TGCT': 4.0, 'C_per_TGCT': 0.0, 
                                        'G_per_TGCT': 0.0, 'T_per_TGCT': 14.0, 'CGCC_per_group': 1.0, 
                                        'CGCC_per_region': 0.6666666666666666, 'CGCT_per_group': 0.0, 
                                        'CGCT_per_region': 0.0, 'TGCC_per_group': 0.0, 'TGCC_per_region': 0.0, 
                                        'TGCT_per_group': 0.5, 'TGCT_per_region': 0.3333333333333333, 
                                        'motif_per_group': 1.5, 'nucleotide_per_group': 9.0, 
                                        'Total_motifs_per_region_length': 12.5, 
                                        'Total_nucleotides_per_region_length': 75.0}]}),
                        # file_stats_dictionary
                        ({"test_files/motif_secondary_features/test1.fasta": 
                                   [{'A_percentage': 22.22222222222222, 'C_percentage': 0.0, 
                                    'G_percentage': 0.0, 'T_percentage': 77.77777777777779, 
                                    'Average_length': 18.0}, 
                                   {'A_percentage': 100.0, 'C_percentage': 0.0, 
                                    'G_percentage': 0.0, 'T_percentage': 0.0, 
                                    'Average_length': 4.0}, 
                                   {'A_percentage': 0.0, 'C_percentage': 0.0, 
                                    'G_percentage': 0.0, 'T_percentage': 100.0, 
                                    'Average_length': 8.0},
                                   {'A_rate_per_motif_per_group': 0.14814814814814814, 
                                    'A_per_group': 2.0, 'C_rate_per_motif_per_group': 0.0, 
                                    'C_per_group': 0.0, 'G_rate_per_motif_per_group': 0.0, 
                                    'G_per_group': 0.0, 'T_rate_per_motif_per_group': 0.5185185185185185, 
                                    'T_per_group': 7.0, 'A_per_CGCC': 2.0, 'C_per_CGCC': 0.0, 
                                    'G_per_CGCC': 0.0, 'T_per_CGCC': 7.0, 'A_per_CGCT': 0.0, 
                                    'C_per_CGCT': 0.0, 'G_per_CGCT': 0.0, 'T_per_CGCT': 0.0, 
                                    'A_per_TGCC': 0.0, 'C_per_TGCC': 0.0, 'G_per_TGCC': 0.0, 
                                    'T_per_TGCC': 0.0, 'A_per_TGCT': 4.0, 'C_per_TGCT': 0.0, 
                                    'G_per_TGCT': 0.0, 'T_per_TGCT': 14.0, 'CGCC_per_group': 1.0, 
                                    'CGCC_per_region': 0.6666666666666666, 'CGCT_per_group': 0.0, 
                                    'CGCT_per_region': 0.0, 'TGCC_per_group': 0.0, 'TGCC_per_region': 0.0, 
                                    'TGCT_per_group': 0.5, 'TGCT_per_region': 0.3333333333333333, 
                                    'motif_per_group': 1.5, 'nucleotide_per_group': 9.0, 
                                    'Total_motifs_per_region_length': 12.5, 
                                    'Total_nucleotides_per_region_length': 75.0}]})),
                        # Test 2: Two identical regions. 
                        # region_dictionary
                        (({"name1": [['AAAATTTTTTTTTTTTTT'], 
                                    ['AAAA'], ['TTTTTTTT'], 
                                    ['TGCTAAAACGCCTTTTTTT', 'TTTTTTTCGCC'], 
                                    24],
                           "name2": [['AAAATTTTTTTTTTTTTT'], 
                                    ['AAAA'], ['TTTTTTTT'], 
                                    ['TGCTAAAACGCCTTTTTTT', 'TTTTTTTCGCC'], 
                                    24]}), 
                        # motif_list
                        (["CGCC", "CGCT", "TGCC", "TGCT"]),
                        # fasta_file_name
                        ("test_files/motif_secondary_features/test2.fasta"),
                        # region_stats_dictionary 
                        {'name1': [{'A': 4, 'C': 0, 'G': 0, 'T': 14, 'Total': 18, 
                                     'A_percentage': 22.22222222222222, 'C_percentage': 0.0, 
                                     'G_percentage': 0.0, 'T_percentage': 77.77777777777779, 
                                     'Average_length': 18.0}, 
                                     {'A': 4, 'C': 0, 'G': 0, 'T': 0, 'Total': 4, 
                                      'A_percentage': 100.0, 'C_percentage': 0.0, 
                                      'G_percentage': 0.0, 'T_percentage': 0.0, 
                                      'Average_length': 4.0}, 
                                      {'A': 0, 'C': 0, 'G': 0, 'T': 8, 'Total': 8, 
                                       'A_percentage': 0.0, 'C_percentage': 0.0, 
                                       'G_percentage': 0.0, 'T_percentage': 100.0, 
                                       'Average_length': 8.0}, 
                                       {'CGCC': 2, 'CGCT': 0, 'TGCC': 0, 'TGCT': 1, 
                                        'Total_motif': 3, 'A': 4, 'C': 0, 'G': 0, 'T': 14, 
                                        'Total_nucleotide': 18, 
                                        'A_rate_per_motif_per_group': 0.14814814814814814, 
                                        'A_per_group': 2.0, 'C_rate_per_motif_per_group': 0.0, 
                                        'C_per_group': 0.0, 'G_rate_per_motif_per_group': 0.0, 
                                        'G_per_group': 0.0, 'T_rate_per_motif_per_group': 0.5185185185185185, 
                                        'T_per_group': 7.0, 'A_per_CGCC': 2.0, 'C_per_CGCC': 0.0, 
                                        'G_per_CGCC': 0.0, 'T_per_CGCC': 7.0, 'A_per_CGCT': 0, 
                                        'C_per_CGCT': 0, 'G_per_CGCT': 0, 'T_per_CGCT': 0, 
                                        'A_per_TGCC': 0, 'C_per_TGCC': 0, 'G_per_TGCC': 0, 
                                        'T_per_TGCC': 0, 'A_per_TGCT': 4.0, 'C_per_TGCT': 0.0, 
                                        'G_per_TGCT': 0.0, 'T_per_TGCT': 14.0, 'CGCC_per_group': 1.0, 
                                        'CGCC_per_region': 0.6666666666666666, 'CGCT_per_group': 0.0, 
                                        'CGCT_per_region': 0.0, 'TGCC_per_group': 0.0, 'TGCC_per_region': 0.0, 
                                        'TGCT_per_group': 0.5, 'TGCT_per_region': 0.3333333333333333, 
                                        'motif_per_group': 1.5, 'nucleotide_per_group': 9.0, 
                                        'Total_motifs_per_region_length': 12.5, 
                                        'Total_nucleotides_per_region_length': 75.0}],
                        'name2': [{'A': 4, 'C': 0, 'G': 0, 'T': 14, 'Total': 18, 
                                     'A_percentage': 22.22222222222222, 'C_percentage': 0.0, 
                                     'G_percentage': 0.0, 'T_percentage': 77.77777777777779, 
                                     'Average_length': 18.0}, 
                                     {'A': 4, 'C': 0, 'G': 0, 'T': 0, 'Total': 4, 
                                      'A_percentage': 100.0, 'C_percentage': 0.0, 
                                      'G_percentage': 0.0, 'T_percentage': 0.0, 
                                      'Average_length': 4.0}, 
                                      {'A': 0, 'C': 0, 'G': 0, 'T': 8, 'Total': 8, 
                                       'A_percentage': 0.0, 'C_percentage': 0.0, 
                                       'G_percentage': 0.0, 'T_percentage': 100.0, 
                                       'Average_length': 8.0}, 
                                       {'CGCC': 2, 'CGCT': 0, 'TGCC': 0, 'TGCT': 1, 
                                        'Total_motif': 3, 'A': 4, 'C': 0, 'G': 0, 'T': 14, 
                                        'Total_nucleotide': 18, 
                                        'A_rate_per_motif_per_group': 0.14814814814814814, 
                                        'A_per_group': 2.0, 'C_rate_per_motif_per_group': 0.0, 
                                        'C_per_group': 0.0, 'G_rate_per_motif_per_group': 0.0, 
                                        'G_per_group': 0.0, 'T_rate_per_motif_per_group': 0.5185185185185185, 
                                        'T_per_group': 7.0, 'A_per_CGCC': 2.0, 'C_per_CGCC': 0.0, 
                                        'G_per_CGCC': 0.0, 'T_per_CGCC': 7.0, 'A_per_CGCT': 0, 
                                        'C_per_CGCT': 0, 'G_per_CGCT': 0, 'T_per_CGCT': 0, 
                                        'A_per_TGCC': 0, 'C_per_TGCC': 0, 'G_per_TGCC': 0, 
                                        'T_per_TGCC': 0, 'A_per_TGCT': 4.0, 'C_per_TGCT': 0.0, 
                                        'G_per_TGCT': 0.0, 'T_per_TGCT': 14.0, 'CGCC_per_group': 1.0, 
                                        'CGCC_per_region': 0.6666666666666666, 'CGCT_per_group': 0.0, 
                                        'CGCT_per_region': 0.0, 'TGCC_per_group': 0.0, 'TGCC_per_region': 0.0, 
                                        'TGCT_per_group': 0.5, 'TGCT_per_region': 0.3333333333333333, 
                                        'motif_per_group': 1.5, 'nucleotide_per_group': 9.0, 
                                        'Total_motifs_per_region_length': 12.5, 
                                        'Total_nucleotides_per_region_length': 75.0}]},
                        # file_stats_dictionary
                        ({"test_files/motif_secondary_features/test2.fasta": 
                                   [{'A_percentage': 22.22222222222222, 'C_percentage': 0.0, 
                                    'G_percentage': 0.0, 'T_percentage': 77.77777777777779, 
                                    'Average_length': 18.0}, 
                                   {'A_percentage': 100.0, 'C_percentage': 0.0, 
                                    'G_percentage': 0.0, 'T_percentage': 0.0, 
                                    'Average_length': 4.0}, 
                                   {'A_percentage': 0.0, 'C_percentage': 0.0, 
                                    'G_percentage': 0.0, 'T_percentage': 100.0, 
                                    'Average_length': 8.0},
                                   {'A_rate_per_motif_per_group': 0.14814814814814814, 
                                    'A_per_group': 2.0, 'C_rate_per_motif_per_group': 0.0, 
                                    'C_per_group': 0.0, 'G_rate_per_motif_per_group': 0.0, 
                                    'G_per_group': 0.0, 'T_rate_per_motif_per_group': 0.5185185185185185, 
                                    'T_per_group': 7.0, 'A_per_CGCC': 2.0, 'C_per_CGCC': 0.0, 
                                    'G_per_CGCC': 0.0, 'T_per_CGCC': 7.0, 'A_per_CGCT': 0.0, 
                                    'C_per_CGCT': 0.0, 'G_per_CGCT': 0.0, 'T_per_CGCT': 0.0, 
                                    'A_per_TGCC': 0.0, 'C_per_TGCC': 0.0, 'G_per_TGCC': 0.0, 
                                    'T_per_TGCC': 0.0, 'A_per_TGCT': 4.0, 'C_per_TGCT': 0.0, 
                                    'G_per_TGCT': 0.0, 'T_per_TGCT': 14.0, 'CGCC_per_group': 1.0, 
                                    'CGCC_per_region': 0.6666666666666666, 'CGCT_per_group': 0.0, 
                                    'CGCT_per_region': 0.0, 'TGCC_per_group': 0.0, 'TGCC_per_region': 0.0, 
                                    'TGCT_per_group': 0.5, 'TGCT_per_region': 0.3333333333333333, 
                                    'motif_per_group': 1.5, 'nucleotide_per_group': 9.0, 
                                    'Total_motifs_per_region_length': 12.5, 
                                    'Total_nucleotides_per_region_length': 75.0}]})),
                        ])

def test_get_motif_counts(region_dictionary, motif_list, fasta_file_name, \
                          region_stats_dictionary, file_stats_dictionary):
    """
        GIVEN a filled region_dictionary that needs to be processed from 
            organize_nucleotides.
        WHEN the function executes other functions to calculate the needed
            statsitics and organize them in two output dictionaries.
        THEN the output dictionaries are checked for the correct data.
    """
    # Executes function and checks output.
    assert m_s_f.get_motif_counts(region_dictionary, motif_list,
                                  fasta_file_name) == (region_stats_dictionary, 
                                                       file_stats_dictionary)

"""
Tests aggregate_motif_counts with various motif_count_list and 
    region size inputs.

    Data Types:
    motif_count_list -- List of dictionaries of motif counts for each
                motif group in a region. Dictionary is organized by motif group.
                [{"CGCC": 0, "CGCT": 0, "TGCC": 0, 
                  "TGCT": 0, "Total_motif": 0,
                  "A": 0, "C": 0, "G": 0, "T": 0, 
                  "Total_nucleotide": 0}].
    region_size -- Integer, size of binding region.
    motif_region_counts -- Dictionary of motif statistics for each region.
        See below for details of entries.
        nucleotide + "rate_per_motif_per_group" -- Nucleotide rate 
        (of total nucleotides) of each nucleotide per motif per group number
        nucleotide + "_per_group" -- Nucleotide type per group.    
        nucleotide + "_per_" + motif -- Nucleotide type per motif type.
        motif + "_per_group" -- Motif fraction type per group.
        motif + "_per_region" -- Motif fraction type per region.
        "motif_per_group" -- Average motifs per group.
        "nucleotide_per_group" -- Average nucleotides per group.
        "Total_motifs_per_region_length" -- Total motifs per region length.
        "Total_nucleotides_per_region" -- Total surrounding nucleotides per region length.
        {A_per_motif_per_group: 0, C_per_motif_per_group: 0, G_per_motif_per_group: 0,
        T_per_motif_per_group: 0, A_per_group: 0, C_per_group: 0, G_per_group: 0,
        T_per_group: 0,
        "A_per_CGCC": 0, "A_per_CGCT": 0, "A_per_TGCC": 0, "A_per_TGCT": 0,
        "C_per_CGCC": 0, "C_per_CGCT": 0, "C_per_TGCC": 0, "C_per_TGCT": 0,
        "G_per_CGCC": 0, "G_per_CGCT": 0, "G_per_TGCC": 0, "G_per_TGCT": 0,
        "T_per_CGCC": 0, "T_per_CGCT": 0, "T_per_TGCC": 0, "T_per_TGCT": 0,
        "CGCC_per_group": 0, "CGCT_per_group": 0, "TGCC_per_group": 0, "TGCT_per_group": 0,
        "CGCC_per_region": 0, "CGCT_per_region": 0, "TGCC_per_region": 0, 
        "TGCT_per_region": 0, "motif_per_group": 0, "nucleotide_per_group": 0,
        "Total_motifs_per_region_length": 0, "Total_nucleotides_per_region_length": 0}.
"""

@pytest.mark.parametrize("motif_count_list, region_size, motif_region_counts", [
                        # Test 1: Simple motif and nucleotide input with one region.
                        # motif_count_list
                        (([{"CGCC": 1, "CGCT": 1, "TGCC": 1,
                            "TGCT": 1, "Total_motif": 4,
                            "A": 2, "C": 1, "G": 1, "T": 6,
                            "Total_nucleotide": 10}]),
                        # region_size
                        (100),
                        # motif_region_counts
                        ({'CGCC': 1, 'CGCT': 1, 'TGCC': 1, 'TGCT': 1, 'Total_motif': 4,
                          'A': 2, 'C': 1, 'G': 1, 'T': 6, 'Total_nucleotide': 10, 
                          'A_rate_per_motif_per_group': 0.05, 'A_per_group': 2.0, 
                          'C_rate_per_motif_per_group': 0.025, 'C_per_group': 1.0, 
                          'G_rate_per_motif_per_group': 0.025, 'G_per_group': 1.0, 
                          'T_rate_per_motif_per_group': 0.15, 'T_per_group': 6.0, 
                          'A_per_CGCC': 2.0, 'C_per_CGCC': 1.0, 'G_per_CGCC': 1.0, 
                          'T_per_CGCC': 6.0, 'A_per_CGCT': 2.0, 'C_per_CGCT': 1.0, 
                          'G_per_CGCT': 1.0, 'T_per_CGCT': 6.0, 'A_per_TGCC': 2.0, 
                          'C_per_TGCC': 1.0, 'G_per_TGCC': 1.0, 'T_per_TGCC': 6.0, 
                          'A_per_TGCT': 2.0, 'C_per_TGCT': 1.0, 'G_per_TGCT': 1.0, 
                          'T_per_TGCT': 6.0, 'CGCC_per_group': 1.0, 'CGCC_per_region': 0.25, 
                          'CGCT_per_group': 1.0, 'CGCT_per_region': 0.25, 'TGCC_per_group': 1.0, 
                          'TGCC_per_region': 0.25, 'TGCT_per_group': 1.0, 'TGCT_per_region': 0.25, 
                          'motif_per_group': 4.0, 'nucleotide_per_group': 10.0, 
                          'Total_motifs_per_region_length': 4.0, 'Total_nucleotides_per_region_length': 10.0})),
                        # Test 2: Simple motif and nucleotide input with three regions.
                        # motif_count_list
                        (([{"CGCC": 1, "CGCT": 1, "TGCC": 1,
                            "TGCT": 1, "Total_motif": 4,
                            "A": 2, "C": 1, "G": 1, "T": 6,
                            "Total_nucleotide": 10},
                            {"CGCC": 1, "CGCT": 1, "TGCC": 1,
                            "TGCT": 1, "Total_motif": 4,
                            "A": 2, "C": 1, "G": 1, "T": 6,
                            "Total_nucleotide": 10},
                            {"CGCC": 1, "CGCT": 1, "TGCC": 1,
                            "TGCT": 1, "Total_motif": 4,
                            "A": 2, "C": 1, "G": 1, "T": 6,
                            "Total_nucleotide": 10}]),
                        # region_size
                        (100),
                        # motif_region_counts
                        ({'CGCC': 3, 'CGCT': 3, 'TGCC': 3, 'TGCT': 3, 'Total_motif': 12,
                          'A': 6, 'C': 3, 'G': 3, 'T': 18, 'Total_nucleotide': 30, 
                          'A_rate_per_motif_per_group': 0.05, 'A_per_group': 2.0, 
                          'C_rate_per_motif_per_group': 0.025, 'C_per_group': 1.0, 
                          'G_rate_per_motif_per_group': 0.025, 'G_per_group': 1.0, 
                          'T_rate_per_motif_per_group': 0.15, 'T_per_group': 6.0, 
                          'A_per_CGCC': 2.0, 'C_per_CGCC': 1.0, 'G_per_CGCC': 1.0, 
                          'T_per_CGCC': 6.0, 'A_per_CGCT': 2.0, 'C_per_CGCT': 1.0, 
                          'G_per_CGCT': 1.0, 'T_per_CGCT': 6.0, 'A_per_TGCC': 2.0, 
                          'C_per_TGCC': 1.0, 'G_per_TGCC': 1.0, 'T_per_TGCC': 6.0, 
                          'A_per_TGCT': 2.0, 'C_per_TGCT': 1.0, 'G_per_TGCT': 1.0, 
                          'T_per_TGCT': 6.0, 'CGCC_per_group': 1.0, 'CGCC_per_region': 0.25, 
                          'CGCT_per_group': 1.0, 'CGCT_per_region': 0.25, 'TGCC_per_group': 1.0, 
                          'TGCC_per_region': 0.25, 'TGCT_per_group': 1.0, 'TGCT_per_region': 0.25, 
                          'motif_per_group': 4.0, 'nucleotide_per_group': 10.0, 
                          'Total_motifs_per_region_length': 12.0, 'Total_nucleotides_per_region_length': 30.0})),
                        # Test 3: Simple motif and nucleotide input with three regions. One
                        # region has zero values and one input has zero values.
                        # motif_count_list
                        (([{"CGCC": 0, "CGCT": 1, "TGCC": 1,
                            "TGCT": 1, "Total_motif": 3,
                            "A": 2, "C": 0, "G": 1, "T": 6,
                            "Total_nucleotide": 9},
                            {"CGCC": 0, "CGCT": 1, "TGCC": 1,
                            "TGCT": 1, "Total_motif": 3,
                            "A": 2, "C": 0, "G": 1, "T": 6,
                            "Total_nucleotide": 9},
                            {"CGCC": 0, "CGCT": 1, "TGCC": 0,
                            "TGCT": 0, "Total_motif": 1,
                            "A": 4, "C": 0, "G": 0, "T": 0,
                            "Total_nucleotide": 4}]),
                        # region_size
                        (100),
                        # motif_region_counts
                        ({'CGCC': 0, 'CGCT': 3, 'TGCC': 2, 'TGCT': 2,
                          'Total_motif': 7, 'A': 8, 'C': 0, 'G': 2, 
                          'T': 12, 'Total_nucleotide': 22, 
                          'A_rate_per_motif_per_group': 0.15584415584415584, 
                          'A_per_group': 2.6666666666666665, 'C_rate_per_motif_per_group': 0.0, 
                          'C_per_group': 0.0, 'G_rate_per_motif_per_group': 0.03896103896103896, 
                          'G_per_group': 0.6666666666666666, 'T_rate_per_motif_per_group': 0.23376623376623373, 
                          'T_per_group': 4.0, 'A_per_CGCC': 0, 'C_per_CGCC': 0, 'G_per_CGCC': 0, 
                          'T_per_CGCC': 0, 'A_per_CGCT': 2.6666666666666665, 'C_per_CGCT': 0.0, 
                          'G_per_CGCT': 0.6666666666666666, 'T_per_CGCT': 4.0, 'A_per_TGCC': 4.0, 
                          'C_per_TGCC': 0.0, 'G_per_TGCC': 1.0, 'T_per_TGCC': 6.0, 'A_per_TGCT': 4.0, 
                          'C_per_TGCT': 0.0, 'G_per_TGCT': 1.0, 'T_per_TGCT': 6.0, 'CGCC_per_group': 0.0, 
                          'CGCC_per_region': 0.0, 'CGCT_per_group': 1.0, 'CGCT_per_region': 0.42857142857142855, 
                          'TGCC_per_group': 0.6666666666666666, 'TGCC_per_region': 0.2857142857142857, 
                          'TGCT_per_group': 0.6666666666666666, 'TGCT_per_region': 0.2857142857142857, 
                          'motif_per_group': 2.3333333333333335, 'nucleotide_per_group': 7.333333333333333, 
                          'Total_motifs_per_region_length': 7.000000000000001, 
                          'Total_nucleotides_per_region_length': 22.0})),     
                        # Test 4: All zero input values.
                        # motif_count_list
                        (([{"CGCC": 0, "CGCT": 0, "TGCC": 0,
                            "TGCT": 0, "Total_motif": 0,
                            "A": 0, "C": 0, "G": 0, "T": 0,
                            "Total_nucleotide": 0}]),
                        # region_size
                        (100),
                        # motif_region_counts
                        ({'CGCC': 0, 'CGCT': 0, 'TGCC': 0, 'TGCT': 0, 'Total_motif': 0,
                          'A': 0, 'C': 0, 'G': 0, 'T': 0, 'Total_nucleotide': 0, 
                          'A_rate_per_motif_per_group': 0, 'A_per_group': 0, 
                          'C_rate_per_motif_per_group': 0, 'C_per_group': 0, 
                          'G_rate_per_motif_per_group': 0, 'G_per_group': 0, 
                          'T_rate_per_motif_per_group': 0, 'T_per_group': 0, 
                          'A_per_CGCC': 0, 'C_per_CGCC': 0, 'G_per_CGCC': 0, 
                          'T_per_CGCC': 0, 'A_per_CGCT': 0, 'C_per_CGCT': 0, 
                          'G_per_CGCT': 0, 'T_per_CGCT': 0, 'A_per_TGCC': 0, 
                          'C_per_TGCC': 0, 'G_per_TGCC': 0, 'T_per_TGCC': 0, 
                          'A_per_TGCT': 0, 'C_per_TGCT': 0, 'G_per_TGCT': 0, 
                          'T_per_TGCT': 0, 'CGCC_per_group': 0, 'CGCC_per_region': 0, 
                          'CGCT_per_group': 0, 'CGCT_per_region': 0, 'TGCC_per_group': 0, 
                          'TGCC_per_region': 0, 'TGCT_per_group': 0, 'TGCT_per_region': 0, 
                          'motif_per_group': 0, 'nucleotide_per_group': 0, 
                          'Total_motifs_per_region_length': 0, 'Total_nucleotides_per_region_length': 0})),
                        ])

def test_aggregate_motif_counts(motif_count_list, region_size, motif_region_counts):
    """
        GIVEN a motif_count_list that has been filled for a binding
            region and the region size.
        WHEN the function totals and averages the motif counts for
            each motif group.
        THEN the motif_region_counts dictionary is checked for the
            correct entries.
    """
    # Executes function and checks output.
    assert m_s_f.aggregate_motif_counts(motif_count_list, region_size) == motif_region_counts

"""
Tests aggregate_file_stats with various motif_count_list and 
    region size inputs.

    Data Types:
    region_stats_dictionary -- Dictionary with motif and seconardary
        statistics for each region. Dictionary is organized by region.
        For nucleotide_count, nucleotide_intra_count, nucleotide_inter_count,
        see nucleotide_content for details, for motif_region_counts, see
        aggregate_motif_counts for details.
        {region: [nucleotide_count, nucleotide_intra_count,
                  nucleotide_inter_count, motif_region_counts]}.    
    file_stats_dictionary -- Dictionary with motif and secondary
        statistics for each input file. Same format as region_stats_dictionary.
        but contains statistics for entire files instead of regions. 
"""

@pytest.mark.parametrize("region_stats_dictionary, file_stats_dictionary", [
                        # Test 1: Single region input.
                        # region_stats_dictionary
                        (({"name1": [{"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.0, "C_percentage": 20.0, 
                                      "G_percentage": 40.0, "T_percentage": 20.0,  
                                      "Total": 10, "Average_length": 10.0},
                                     {"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.0, "C_percentage": 20.0, 
                                      "G_percentage": 40.0, "T_percentage": 20.0,  
                                      "Total": 10, "Average_length": 10.0},
                                     {"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.0, "C_percentage": 20.0, 
                                      "G_percentage": 40.0, "T_percentage": 20.0,  
                                      "Total": 10, "Average_length": 10.0},
                                     {'CGCC': 1, 'CGCT': 1, 'TGCC': 1, 'TGCT': 1, 'Total_motif': 4,
                                      'A': 2, 'C': 1, 'G': 1, 'T': 6, 'Total_nucleotide': 10, 
                                      'A_rate_per_motif_per_group': 0.05, 'A_per_group': 2.0, 
                                      'C_rate_per_motif_per_group': 0.025, 'C_per_group': 1.0, 
                                      'G_rate_per_motif_per_group': 0.025, 'G_per_group': 1.0, 
                                      'T_rate_per_motif_per_group': 0.15, 'T_per_group': 6.0, 
                                      'A_per_CGCC': 2.0, 'C_per_CGCC': 1.0, 'G_per_CGCC': 1.0, 
                                      'T_per_CGCC': 6.0, 'A_per_CGCT': 2.0, 'C_per_CGCT': 1.0, 
                                      'G_per_CGCT': 1.0, 'T_per_CGCT': 6.0, 'A_per_TGCC': 2.0, 
                                      'C_per_TGCC': 1.0, 'G_per_TGCC': 1.0, 'T_per_TGCC': 6.0, 
                                      'A_per_TGCT': 2.0, 'C_per_TGCT': 1.0, 'G_per_TGCT': 1.0, 
                                      'T_per_TGCT': 6.0, 'CGCC_per_group': 1.0, 'CGCC_per_region': 0.25, 
                                      'CGCT_per_group': 1.0, 'CGCT_per_region': 0.25, 'TGCC_per_group': 1.0, 
                                      'TGCC_per_region': 0.25, 'TGCT_per_group': 1.0, 'TGCT_per_region': 0.25, 
                                      'motif_per_group': 4.0, 'nucleotide_per_group': 10.0, 
                                      'Total_motifs_per_region_length': 4.0, 'Total_nucleotides_per_region_length': 10.0}]}),
                          # file_stats_dictionary
                          ({"file": [{"A_percentage": 20.0, "C_percentage": 20.0, 
                                      "G_percentage": 40.0, "T_percentage": 20.0,  
                                      "Average_length": 10.0},
                                     {"A_percentage": 20.0, "C_percentage": 20.0, 
                                      "G_percentage": 40.0, "T_percentage": 20.0,  
                                      "Average_length": 10.0},
                                     {"A_percentage": 20.0, "C_percentage": 20.0, 
                                      "G_percentage": 40.0, "T_percentage": 20.0,  
                                      "Average_length": 10.0},
                                     {'A_rate_per_motif_per_group': 0.05, 'A_per_group': 2.0, 
                                      'C_rate_per_motif_per_group': 0.025, 'C_per_group': 1.0, 
                                      'G_rate_per_motif_per_group': 0.025, 'G_per_group': 1.0, 
                                      'T_rate_per_motif_per_group': 0.15, 'T_per_group': 6.0, 
                                      'A_per_CGCC': 2.0, 'C_per_CGCC': 1.0, 'G_per_CGCC': 1.0, 
                                      'T_per_CGCC': 6.0, 'A_per_CGCT': 2.0, 'C_per_CGCT': 1.0, 
                                      'G_per_CGCT': 1.0, 'T_per_CGCT': 6.0, 'A_per_TGCC': 2.0, 
                                      'C_per_TGCC': 1.0, 'G_per_TGCC': 1.0, 'T_per_TGCC': 6.0, 
                                      'A_per_TGCT': 2.0, 'C_per_TGCT': 1.0, 'G_per_TGCT': 1.0, 
                                      'T_per_TGCT': 6.0, 'CGCC_per_group': 1.0, 'CGCC_per_region': 0.25, 
                                      'CGCT_per_group': 1.0, 'CGCT_per_region': 0.25, 'TGCC_per_group': 1.0, 
                                      'TGCC_per_region': 0.25, 'TGCT_per_group': 1.0, 'TGCT_per_region': 0.25, 
                                      'motif_per_group': 4.0, 'nucleotide_per_group': 10.0, 
                                      'Total_motifs_per_region_length': 4.0, 'Total_nucleotides_per_region_length': 10.0}]})),
                          # Test 2: Three region input.
                          # region_stats_dictionary
                          (({"name1": [{"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.0, "C_percentage": 20.0, 
                                      "G_percentage": 40.0, "T_percentage": 20.0,  
                                      "Total": 10, "Average_length": 10.0},  
                                     {"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.0, "C_percentage": 20.0, 
                                      "G_percentage": 40.0, "T_percentage": 20.0,  
                                      "Total": 10, "Average_length": 10.0},
                                     {"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.0, "C_percentage": 20.0, 
                                      "G_percentage": 40.0, "T_percentage": 20.0,  
                                      "Total": 10, "Average_length": 10.0},
                                     {'CGCC': 1, 'CGCT': 1, 'TGCC': 1, 'TGCT': 1, 'Total_motif': 4,
                                      'A': 2, 'C': 1, 'G': 1, 'T': 6, 'Total_nucleotide': 10, 
                                      'A_rate_per_motif_per_group': 1.0, 'A_per_group': 1.0, 
                                      'C_rate_per_motif_per_group': 1.0, 'C_per_group': 1.0, 
                                      'G_rate_per_motif_per_group': 1.0, 'G_per_group': 1.0, 
                                      'T_rate_per_motif_per_group': 1.0, 'T_per_group': 1.0, 
                                      'A_per_CGCC': 1.0, 'C_per_CGCC': 1.0, 'G_per_CGCC': 1.0, 
                                      'T_per_CGCC': 1.0, 'A_per_CGCT': 1.0, 'C_per_CGCT': 1.0, 
                                      'G_per_CGCT': 1.0, 'T_per_CGCT': 1.0, 'A_per_TGCC': 1.0, 
                                      'C_per_TGCC': 1.0, 'G_per_TGCC': 1.0, 'T_per_TGCC': 1.0, 
                                      'A_per_TGCT': 1.0, 'C_per_TGCT': 1.0, 'G_per_TGCT': 1.0, 
                                      'T_per_TGCT': 1.0, 'CGCC_per_group': 1.0, 'CGCC_per_region': 1.0, 
                                      'CGCT_per_group': 1.0, 'CGCT_per_region': 1.0, 'TGCC_per_group': 1.0, 
                                      'TGCC_per_region': 1.0, 'TGCT_per_group': 1.0, 'TGCT_per_region': 1.0, 
                                      'motif_per_group': 1.0, 'nucleotide_per_group': 1.0, 
                                      'Total_motifs_per_region_length': 1.0, 'Total_nucleotides_per_region_length': 1.0}],
                            "name2": [{"A": 4, "C": 4, "G": 6, "T": 4,
                                      "A_percentage": 40.0, "C_percentage": 40.0, 
                                      "G_percentage": 60.0, "T_percentage": 40.0,  
                                      "Total": 10, "Average_length": 20.0},  
                                     {"A": 4, "C": 4, "G": 6, "T": 4,
                                      "A_percentage": 40.0, "C_percentage": 40.0, 
                                      "G_percentage": 60.0, "T_percentage": 40.0,  
                                      "Total": 10, "Average_length": 20.0},
                                     {"A": 4, "C": 4, "G": 6, "T": 4,
                                      "A_percentage": 40.0, "C_percentage": 40.0, 
                                      "G_percentage": 60.0, "T_percentage": 40.0,  
                                      "Total": 10, "Average_length": 20.0},
                                      {'CGCC': 2, 'CGCT': 2, 'TGCC': 2, 'TGCT': 2, 'Total_motif': 8,
                                      'A': 4, 'C': 2, 'G': 4, 'T': 8, 'Total_nucleotide': 18, 
                                      'A_rate_per_motif_per_group': 2.0, 'A_per_group': 2.0, 
                                      'C_rate_per_motif_per_group': 2.0, 'C_per_group': 2.0, 
                                      'G_rate_per_motif_per_group': 2.0, 'G_per_group': 2.0, 
                                      'T_rate_per_motif_per_group': 2.0, 'T_per_group': 2.0, 
                                      'A_per_CGCC': 2.0, 'C_per_CGCC': 2.0, 'G_per_CGCC': 2.0, 
                                      'T_per_CGCC': 2.0, 'A_per_CGCT': 2.0, 'C_per_CGCT': 2.0, 
                                      'G_per_CGCT': 2.0, 'T_per_CGCT': 2.0, 'A_per_TGCC': 2.0, 
                                      'C_per_TGCC': 2.0, 'G_per_TGCC': 2.0, 'T_per_TGCC': 2.0, 
                                      'A_per_TGCT': 2.0, 'C_per_TGCT': 2.0, 'G_per_TGCT': 2.0, 
                                      'T_per_TGCT': 2.0, 'CGCC_per_group': 2.0, 'CGCC_per_region': 2.0, 
                                      'CGCT_per_group': 2.0, 'CGCT_per_region': 2.0, 'TGCC_per_group': 2.0, 
                                      'TGCC_per_region': 2.0, 'TGCT_per_group': 2.0, 'TGCT_per_region': 2.0, 
                                      'motif_per_group': 2.0, 'nucleotide_per_group': 2.0, 
                                      'Total_motifs_per_region_length': 2.0, 'Total_nucleotides_per_region_length': 2.0}],
                            "name3": [{"A": 6, "C": 6, "G": 8, "T": 6,
                                      "A_percentage": 60.0, "C_percentage": 60.0, 
                                      "G_percentage": 80.0, "T_percentage": 60.0,  
                                      "Total": 10, "Average_length": 30.0},  
                                     {"A": 6, "C": 6, "G": 8, "T": 6,
                                      "A_percentage": 60.0, "C_percentage": 60.0, 
                                      "G_percentage": 80.0, "T_percentage": 60.0,  
                                      "Total": 10, "Average_length": 30.0},
                                     {"A": 6, "C": 6, "G": 8, "T": 6,
                                      "A_percentage": 60.0, "C_percentage": 60.0, 
                                      "G_percentage": 80.0, "T_percentage": 60.0,  
                                      "Total": 10, "Average_length": 30.0},
                                      {'CGCC': 3, 'CGCT': 3, 'TGCC': 3, 'TGCT': 3, 'Total_motif': 12,
                                      'A': 6, 'C': 3, 'G': 6, 'T': 10, 'Total_nucleotide': 25, 
                                      'A_rate_per_motif_per_group': 3.0, 'A_per_group': 3.0, 
                                      'C_rate_per_motif_per_group': 3.0, 'C_per_group': 3.0, 
                                      'G_rate_per_motif_per_group': 3.0, 'G_per_group': 3.0, 
                                      'T_rate_per_motif_per_group': 3.0, 'T_per_group': 3.0, 
                                      'A_per_CGCC': 3.0, 'C_per_CGCC': 3.0, 'G_per_CGCC': 3.0, 
                                      'T_per_CGCC': 3.0, 'A_per_CGCT': 3.0, 'C_per_CGCT': 3.0, 
                                      'G_per_CGCT': 3.0, 'T_per_CGCT': 3.0, 'A_per_TGCC': 3.0, 
                                      'C_per_TGCC': 3.0, 'G_per_TGCC': 3.0, 'T_per_TGCC': 3.0, 
                                      'A_per_TGCT': 3.0, 'C_per_TGCT': 3.0, 'G_per_TGCT': 3.0, 
                                      'T_per_TGCT': 3.0, 'CGCC_per_group': 3.0, 'CGCC_per_region': 3.0, 
                                      'CGCT_per_group': 3.0, 'CGCT_per_region': 3.0, 'TGCC_per_group': 3.0, 
                                      'TGCC_per_region': 3.0, 'TGCT_per_group': 3.0, 'TGCT_per_region': 3.0, 
                                      'motif_per_group': 3.0, 'nucleotide_per_group': 3.0, 
                                      'Total_motifs_per_region_length': 3.0, 'Total_nucleotides_per_region_length': 3.0}]}),
                        # file_stats_dictionary
                        ({"file": [{"A_percentage": 40.0, "C_percentage": 40.0, 
                                      "G_percentage": 60.0, "T_percentage": 40.0,  
                                      "Average_length": 20.0},
                                     {"A_percentage": 40.0, "C_percentage": 40.0, 
                                      "G_percentage": 60.0, "T_percentage": 40.0,  
                                      "Average_length": 20.0},
                                     {"A_percentage": 40.0, "C_percentage": 40.0, 
                                      "G_percentage": 60.0, "T_percentage": 40.0,  
                                      "Average_length": 20.0},
                                     {'A_rate_per_motif_per_group': 2.0, 'A_per_group': 2.0, 
                                      'C_rate_per_motif_per_group': 2.0, 'C_per_group': 2.0, 
                                      'G_rate_per_motif_per_group': 2.0, 'G_per_group': 2.0, 
                                      'T_rate_per_motif_per_group': 2.0, 'T_per_group': 2.0, 
                                      'A_per_CGCC': 2.0, 'C_per_CGCC': 2.0, 'G_per_CGCC': 2.0, 
                                      'T_per_CGCC': 2.0, 'A_per_CGCT': 2.0, 'C_per_CGCT': 2.0, 
                                      'G_per_CGCT': 2.0, 'T_per_CGCT': 2.0, 'A_per_TGCC': 2.0, 
                                      'C_per_TGCC': 2.0, 'G_per_TGCC': 2.0, 'T_per_TGCC': 2.0, 
                                      'A_per_TGCT': 2.0, 'C_per_TGCT': 2.0, 'G_per_TGCT': 2.0, 
                                      'T_per_TGCT': 2.0, 'CGCC_per_group': 2.0, 'CGCC_per_region': 2.0, 
                                      'CGCT_per_group': 2.0, 'CGCT_per_region': 2.0, 'TGCC_per_group': 2.0, 
                                      'TGCC_per_region': 2.0, 'TGCT_per_group': 2.0, 'TGCT_per_region': 2.0, 
                                      'motif_per_group': 2.0, 'nucleotide_per_group': 2.0, 
                                      'Total_motifs_per_region_length': 2.0, 'Total_nucleotides_per_region_length': 2.0}]})),
                        # Test 3: Zero values.
                        # region_stats_dictionary
                        (({"name1": [{"A": 0, "C": 0, "G": 0, "T": 0,
                                      "A_percentage": 0, "C_percentage": 0, 
                                      "G_percentage": 0, "T_percentage": 0,  
                                      "Total": 0, "Average_length": 0},
                                     {"A": 0, "C": 0, "G": 0, "T": 0,
                                      "A_percentage": 0, "C_percentage": 0, 
                                      "G_percentage": 0, "T_percentage": 0,  
                                      "Total": 0, "Average_length": 0},
                                     {"A": 0, "C": 0, "G": 0, "T": 0,
                                      "A_percentage": 0, "C_percentage": 0, 
                                      "G_percentage": 0, "T_percentage": 0,  
                                      "Total": 0, "Average_length": 0},
                                     {'CGCC': 0, 'CGCT': 0, 'TGCC': 0, 'TGCT': 0, 'Total_motif': 0,
                                      'A': 0, 'C': 0, 'G': 0, 'T': 0, 'Total_nucleotide': 0, 
                                      'A_rate_per_motif_per_group': 0, 'A_per_group': 0, 
                                      'C_rate_per_motif_per_group': 0, 'C_per_group': 0, 
                                      'G_rate_per_motif_per_group': 0, 'G_per_group': 0, 
                                      'T_rate_per_motif_per_group': 0, 'T_per_group': 0, 
                                      'A_per_CGCC': 0, 'C_per_CGCC': 0, 'G_per_CGCC': 0, 
                                      'T_per_CGCC': 0, 'A_per_CGCT': 0, 'C_per_CGCT': 0, 
                                      'G_per_CGCT': 0, 'T_per_CGCT': 0, 'A_per_TGCC': 0, 
                                      'C_per_TGCC': 0, 'G_per_TGCC': 0, 'T_per_TGCC': 0, 
                                      'A_per_TGCT': 0, 'C_per_TGCT': 0, 'G_per_TGCT': 0, 
                                      'T_per_TGCT': 0, 'CGCC_per_group': 0, 'CGCC_per_region': 0, 
                                      'CGCT_per_group': 0, 'CGCT_per_region': 0, 'TGCC_per_group': 0, 
                                      'TGCC_per_region': 0, 'TGCT_per_group': 0, 'TGCT_per_region': 0, 
                                      'motif_per_group': 0, 'nucleotide_per_group': 0, 
                                      'Total_motifs_per_region_length': 0, 'Total_nucleotides_per_region_length': 0}]}),
                        # file_stats_dictionary
                        ({"file": [{"A_percentage": 0.0, "C_percentage": 0.0, 
                                      "G_percentage": 0.0, "T_percentage": 0.0,  
                                      "Average_length": 0.0},
                                     {"A_percentage": 0.0, "C_percentage": 0.0, 
                                      "G_percentage": 0.0, "T_percentage": 0.0,  
                                      "Average_length": 0.0},
                                     {"A_percentage": 0.0, "C_percentage": 0.0, 
                                      "G_percentage": 0.0, "T_percentage": 0.0,  
                                      "Average_length": 0.0},
                                     {'A_rate_per_motif_per_group': 0.0, 'A_per_group': 0.0, 
                                      'C_rate_per_motif_per_group': 0.0, 'C_per_group': 0.0, 
                                      'G_rate_per_motif_per_group': 0.0, 'G_per_group': 0.0, 
                                      'T_rate_per_motif_per_group': 0.0, 'T_per_group': 0.0, 
                                      'A_per_CGCC': 0.0, 'C_per_CGCC': 0.0, 'G_per_CGCC': 0.0, 
                                      'T_per_CGCC': 0.0, 'A_per_CGCT': 0.0, 'C_per_CGCT': 0.0, 
                                      'G_per_CGCT': 0.0, 'T_per_CGCT': 0.0, 'A_per_TGCC': 0.0, 
                                      'C_per_TGCC': 0.0, 'G_per_TGCC': 0.0, 'T_per_TGCC': 0.0, 
                                      'A_per_TGCT': 0.0, 'C_per_TGCT': 0.0, 'G_per_TGCT': 0.0, 
                                      'T_per_TGCT': 0.0, 'CGCC_per_group': 0.0, 'CGCC_per_region': 0.0, 
                                      'CGCT_per_group': 0.0, 'CGCT_per_region': 0.0, 'TGCC_per_group': 0.0, 
                                      'TGCC_per_region': 0.0, 'TGCT_per_group': 0.0, 'TGCT_per_region': 0.0, 
                                      'motif_per_group': 0.0, 'nucleotide_per_group': 0.0, 
                                      'Total_motifs_per_region_length': 0.0, 'Total_nucleotides_per_region_length': 0.0}]})),
                             #Test 4: Empty input.
                             # region_stats_dictionary
                             ({},
                              # file_stats_dictionary
                              {'file': [{}, {}, {}, {}]}), 
                                ])

def test_aggregate_file_stats(region_stats_dictionary, file_stats_dictionary):
    """
        GIVEN a region_stats_dictionary that was generated by get_motif_counts
            that needs file counts.
        WHEN the function totals and averages each of the statisitics in the
            dictionary on a per file basis and adds to a dictionary.
        THEN the output dictionary is checked for the correct data.
    """
    # Executes function and checks output.
    assert m_s_f.aggregate_file_stats(region_stats_dictionary, "file") == file_stats_dictionary

"""
Tests motif_content with various motif_count_list and 
    region size inputs.

    Data Types:
    motifs_group_list -- List of motifs and surrounding nucleotides.
    motif_list -- List of motifs to be analyzed. Default is YGCY.
    region_size -- Integer, size of binding region.
    motif_region_counts -- Dictionary of motif statistics for each region.
        See below for details of entries.
        nucleotide + "rate_per_motif_per_group" -- Nucleotide rate 
        (of total nucleotides) of each nucleotide per motif per group number
        nucleotide + "_per_group" -- Nucleotide type per group.    
        nucleotide + "_per_" + motif -- Nucleotide type per motif type.
        motif + "_per_group" -- Motif fraction type per group.
        motif + "_per_region" -- Motif fraction type per region.
        "motif_per_group" -- Average motifs per group.
        "nucleotide_per_group" -- Average nucleotides per group.
        "Total_motifs_per_region_length" -- Total motifs per region length.
        "Total_nucleotides_per_region" -- Total surrounding nucleotides per region length.
        {A_per_motif_per_group: 0, C_per_motif_per_group: 0, G_per_motif_per_group: 0,
        T_per_motif_per_group: 0, A_per_group: 0, C_per_group: 0, G_per_group: 0,
        T_per_group: 0,
        "A_per_CGCC": 0, "A_per_CGCT": 0, "A_per_TGCC": 0, "A_per_TGCT": 0,
        "C_per_CGCC": 0, "C_per_CGCT": 0, "C_per_TGCC": 0, "C_per_TGCT": 0,
        "G_per_CGCC": 0, "G_per_CGCT": 0, "G_per_TGCC": 0, "G_per_TGCT": 0,
        "T_per_CGCC": 0, "T_per_CGCT": 0, "T_per_TGCC": 0, "T_per_TGCT": 0,
        "CGCC_per_group": 0, "CGCT_per_group": 0, "TGCC_per_group": 0, "TGCT_per_group": 0,
        "CGCC_per_region": 0, "CGCT_per_region": 0, "TGCC_per_region": 0, 
        "TGCT_per_region": 0, "motif_per_group": 0, "nucleotide_per_group": 0,
        "Total_motifs_per_region_length": 0, "Total_nucleotides_per_region_length": 0}.
"""

@pytest.mark.parametrize("motifs_group_list, motif_list, region_size, motif_region_counts", [
                        # Test 1: Simple motif and nucleotide input with one region.
                        # motifs_group_list
                        ((["CACGCCTTCGCTTTTGCCTTTGCTAG"]),
                        # motif_list
                        (["CGCC", "CGCT", "TGCC", "TGCT"]),  
                        # region_size
                        (100),
                        # motif_region_counts
                        ({'CGCC': 1, 'CGCT': 1, 'TGCC': 1, 'TGCT': 1, 'Total_motif': 4,
                          'A': 2, 'C': 1, 'G': 1, 'T': 6, 'Total_nucleotide': 10, 
                          'A_rate_per_motif_per_group': 0.05, 'A_per_group': 2.0, 
                          'C_rate_per_motif_per_group': 0.025, 'C_per_group': 1.0, 
                          'G_rate_per_motif_per_group': 0.025, 'G_per_group': 1.0, 
                          'T_rate_per_motif_per_group': 0.15, 'T_per_group': 6.0, 
                          'A_per_CGCC': 2.0, 'C_per_CGCC': 1.0, 'G_per_CGCC': 1.0, 
                          'T_per_CGCC': 6.0, 'A_per_CGCT': 2.0, 'C_per_CGCT': 1.0, 
                          'G_per_CGCT': 1.0, 'T_per_CGCT': 6.0, 'A_per_TGCC': 2.0, 
                          'C_per_TGCC': 1.0, 'G_per_TGCC': 1.0, 'T_per_TGCC': 6.0, 
                          'A_per_TGCT': 2.0, 'C_per_TGCT': 1.0, 'G_per_TGCT': 1.0, 
                          'T_per_TGCT': 6.0, 'CGCC_per_group': 1.0, 'CGCC_per_region': 0.25, 
                          'CGCT_per_group': 1.0, 'CGCT_per_region': 0.25, 'TGCC_per_group': 1.0, 
                          'TGCC_per_region': 0.25, 'TGCT_per_group': 1.0, 'TGCT_per_region': 0.25, 
                          'motif_per_group': 4.0, 'nucleotide_per_group': 10.0, 
                          'Total_motifs_per_region_length': 4.0, 'Total_nucleotides_per_region_length': 10.0})),
                        # Test 2: Simple motif and nucleotide input with three regions.
                        # motifs_group_list
                        ((["CACGCCTTCGCTTTTGCCTTTGCTAG", "CACGCCTTCGCTTTTGCCTTTGCTAG",
                            "CACGCCTTCGCTTTTGCCTTTGCTAG"]),
                        # motif_list
                        (["CGCC", "CGCT", "TGCC", "TGCT"]),    
                        # region_size
                        (100),
                        # motif_region_counts
                        ({'CGCC': 3, 'CGCT': 3, 'TGCC': 3, 'TGCT': 3, 'Total_motif': 12,
                          'A': 6, 'C': 3, 'G': 3, 'T': 18, 'Total_nucleotide': 30, 
                          'A_rate_per_motif_per_group': 0.05, 'A_per_group': 2.0, 
                          'C_rate_per_motif_per_group': 0.025, 'C_per_group': 1.0, 
                          'G_rate_per_motif_per_group': 0.025, 'G_per_group': 1.0, 
                          'T_rate_per_motif_per_group': 0.15, 'T_per_group': 6.0, 
                          'A_per_CGCC': 2.0, 'C_per_CGCC': 1.0, 'G_per_CGCC': 1.0, 
                          'T_per_CGCC': 6.0, 'A_per_CGCT': 2.0, 'C_per_CGCT': 1.0, 
                          'G_per_CGCT': 1.0, 'T_per_CGCT': 6.0, 'A_per_TGCC': 2.0, 
                          'C_per_TGCC': 1.0, 'G_per_TGCC': 1.0, 'T_per_TGCC': 6.0, 
                          'A_per_TGCT': 2.0, 'C_per_TGCT': 1.0, 'G_per_TGCT': 1.0, 
                          'T_per_TGCT': 6.0, 'CGCC_per_group': 1.0, 'CGCC_per_region': 0.25, 
                          'CGCT_per_group': 1.0, 'CGCT_per_region': 0.25, 'TGCC_per_group': 1.0, 
                          'TGCC_per_region': 0.25, 'TGCT_per_group': 1.0, 'TGCT_per_region': 0.25, 
                          'motif_per_group': 4.0, 'nucleotide_per_group': 10.0, 
                          'Total_motifs_per_region_length': 12.0, 'Total_nucleotides_per_region_length': 30.0})),
                        # Test 3: Simple motif and nucleotide input with three regions. One
                        # region has zero values and one nucleotide has zero values. Also
                        # tests for multiple motifs in a row and motifs that come up in the
                        # first part of the sequence and the last part of the sequence.
                        # motifs_group_list
                        ((["TTTCGCTTGCCTGCTAGATTT", "CGCTTGCCTTTATTTTGCTAG",
                           "AAAACGCT"]), 
                        # motif_list
                        (["CGCC", "CGCT", "TGCC", "TGCT"]),    
                        # region_size
                        (100),
                        # motif_region_counts
                        ({'CGCC': 0, 'CGCT': 3, 'TGCC': 2, 'TGCT': 2,
                          'Total_motif': 7, 'A': 8, 'C': 0, 'G': 2, 
                          'T': 12, 'Total_nucleotide': 22, 
                          'A_rate_per_motif_per_group': 0.15584415584415584, 
                          'A_per_group': 2.6666666666666665, 'C_rate_per_motif_per_group': 0.0, 
                          'C_per_group': 0.0, 'G_rate_per_motif_per_group': 0.03896103896103896, 
                          'G_per_group': 0.6666666666666666, 'T_rate_per_motif_per_group': 0.23376623376623373, 
                          'T_per_group': 4.0, 'A_per_CGCC': 0, 'C_per_CGCC': 0, 'G_per_CGCC': 0, 
                          'T_per_CGCC': 0, 'A_per_CGCT': 2.6666666666666665, 'C_per_CGCT': 0.0, 
                          'G_per_CGCT': 0.6666666666666666, 'T_per_CGCT': 4.0, 'A_per_TGCC': 4.0, 
                          'C_per_TGCC': 0.0, 'G_per_TGCC': 1.0, 'T_per_TGCC': 6.0, 'A_per_TGCT': 4.0, 
                          'C_per_TGCT': 0.0, 'G_per_TGCT': 1.0, 'T_per_TGCT': 6.0, 'CGCC_per_group': 0.0, 
                          'CGCC_per_region': 0.0, 'CGCT_per_group': 1.0, 'CGCT_per_region': 0.42857142857142855, 
                          'TGCC_per_group': 0.6666666666666666, 'TGCC_per_region': 0.2857142857142857, 
                          'TGCT_per_group': 0.6666666666666666, 'TGCT_per_region': 0.2857142857142857, 
                          'motif_per_group': 2.3333333333333335, 'nucleotide_per_group': 7.333333333333333, 
                          'Total_motifs_per_region_length': 7.000000000000001, 
                          'Total_nucleotides_per_region_length': 22.0})),     
                        # Test 4: All zero input values.
                        # motifs_group_list
                        (([]),
                        # motif_list
                        (["CGCC", "CGCT", "TGCC", "TGCT"]),    
                        # region_size
                        (100),
                        # motif_region_counts
                        ({'CGCC': 0, 'CGCT': 0, 'TGCC': 0, 'TGCT': 0, 'Total_motif': 0,
                          'A': 0, 'C': 0, 'G': 0, 'T': 0, 'Total_nucleotide': 0, 
                          'A_rate_per_motif_per_group': 0, 'A_per_group': 0, 
                          'C_rate_per_motif_per_group': 0, 'C_per_group': 0, 
                          'G_rate_per_motif_per_group': 0, 'G_per_group': 0, 
                          'T_rate_per_motif_per_group': 0, 'T_per_group': 0, 
                          'A_per_CGCC': 0, 'C_per_CGCC': 0, 'G_per_CGCC': 0, 
                          'T_per_CGCC': 0, 'A_per_CGCT': 0, 'C_per_CGCT': 0, 
                          'G_per_CGCT': 0, 'T_per_CGCT': 0, 'A_per_TGCC': 0, 
                          'C_per_TGCC': 0, 'G_per_TGCC': 0, 'T_per_TGCC': 0, 
                          'A_per_TGCT': 0, 'C_per_TGCT': 0, 'G_per_TGCT': 0, 
                          'T_per_TGCT': 0, 'CGCC_per_group': 0, 'CGCC_per_region': 0, 
                          'CGCT_per_group': 0, 'CGCT_per_region': 0, 'TGCC_per_group': 0, 
                          'TGCC_per_region': 0, 'TGCT_per_group': 0, 'TGCT_per_region': 0, 
                          'motif_per_group': 0, 'nucleotide_per_group': 0, 
                          'Total_motifs_per_region_length': 0, 'Total_nucleotides_per_region_length': 0})),
                        ])

def test_motif_content(motifs_group_list, motif_list, 
                       region_size, motif_region_counts):
    """
        GIVEN a list of motifs that needs to be identified (motif_list), 
            a motifs_group_list, and the region size that needs to 
            be analyzed for motif content.
        WHEN the function totals and averages the motif counts for
            each motif group.
        THEN the motif_region_counts dictionary is checked for the
            correct entries.
    """
    # Executes function and checks output.
    assert m_s_f.motif_content(motifs_group_list, motif_list, 
                               region_size) == motif_region_counts

"""
Tests nucleotide_content with various nucleotide lists.

    Data Types:
    nucleotide_list -- List of nucleotide strings to be analyzed.
    nucleotide_count -- Dictionary of nucleotide counts for each nucleotide
        in nucleotide_list. Dictionary is organized by nucleotide. Contains
        nucleotide count total, nucleotide percentage, and average nucleotide
        length.
        nucleotide_count = {"A": 2, "C": 2, "G": 4, "T": 2,
                          "A_percentage": 0.0, "C_percentage": 0.0, 
                          "G_percentage": 0.0, "T_percentage": 0.0,  
                          "Total": 0, "Average_length": 0.0}.
"""

@pytest.mark.parametrize("nucleotide_list, nucleotide_count", [
                        # Test 1: Single single nucleotide string input.
                        # nucleotide_list
                        ((["AAGGTTCCGG"]),
                        # nucleotide_count
                        ({"A": 2, "C": 2, "G": 4, "T": 2,
                          "A_percentage": 20.0, "C_percentage": 20.0, 
                          "G_percentage": 40.0, "T_percentage": 20.0,  
                          "Total": 10, "Average_length": 10.0})),
                        # Test 2: Multiple single nucleotide string input.
                        # nucleotide_list
                        ((["AAGGTTCCGG", "AAGGTTCCGG", "AATTCC"]),
                        # nucleotide_count
                        ({"A": 6, "C": 6, "G": 8, "T": 6,
                          "A_percentage": 23.076923076923077, "C_percentage": 23.076923076923077, 
                          "G_percentage": 30.76923076923077, "T_percentage": 23.076923076923077,  
                          "Total": 26, "Average_length": 8.666666666666666})),
                        # Test 3: Multiple single nucleotide string input, one with a single 
                        # nucleotide and one string that has zero nucleotides.
                        # nucleotide_list
                        ((["AAGGTTCCGG", "AAGGTTCCGG", "AATTCC", "", "A"]),
                        # nucleotide_count
                        ({"A": 7, "C": 6, "G": 8, "T": 6,
                          "A_percentage": 25.925925925925924, "C_percentage": 22.22222222222222, 
                          "G_percentage": 29.629629629629626, "T_percentage": 22.22222222222222,  
                          "Total": 27,"Average_length": 5.4})),
                        # Test 4: Empty list.
                        # nucleotide_list
                        (([]),
                        # nucleotide_count
                        ({"A": 0, "C": 0, "G": 0, "T": 0,
                          "A_percentage": 0.0, "C_percentage": 0.0, 
                          "G_percentage": 0.0, "T_percentage": 0.0,  
                          "Total": 0.0,"Average_length": 0.0})),                 
                        ])

def test_nucleotide_content(nucleotide_list, nucleotide_count):
    """
        GIVEN a list of nucleotide strings that need to be analyzed.
        WHEN the function totals and averages the nucleotide counts
            for all strings in the list.
        THEN the nucleotide_count dictionary is checked for the
            correct entries.
    """
    # Executes function and checks output.
    assert m_s_f.nucleotide_content(nucleotide_list) == nucleotide_count

"""
Tests organize_files with various FASTA files. Kmer size and different
    motifs not tested.

    Data Types:
    region_fasta_list -- List of FASTA files of RBP binding regions.
        should include previously generated random regions from
        make_random_bed_for_motif.py. The random files are detected
        by searching for "random" in the file name.
    kmer_size -- Size of kmer to be analyzed.
    motif_list -- List of motifs to be analyzed. Default is YGCY.
    motif_title -- Title for motif output file. Default is YGCY.
    output_file_prefix -- Prefix for output files.
    region_output_check_list -- List of CSV files of expected 
        nucleotide and motif secondary statistics.
    count_match_list  -- List of number of lines that match 
        between the output and the first check file. 
"""

@pytest.mark.parametrize("region_fasta_list, kmer_size, motif_list, motif_title, \
                          output_file_prefix, region_output_check_list, count_match_list", [
                        # Test 1: Single region with two motif groups. One
                        # group with two motifs, the other with one.
                        # region_fasta_list
                        ((["test_files/motif_secondary_features/test1.fasta"]),
                        # kmer_size
                        (4),
                        # motif_list
                        (["CGCC", "CGCT", "TGCC", "TGCT"]),
                        # motif_title
                        ("YGCY"),
                        # output_file_prefix
                        ("test_files/motif_secondary_features/test12"),
                        # region_output_check_list
                        (["test_files/motif_secondary_features/test12_secondary_motif_stats_per_region_1_check.csv",
                          "test_files/motif_secondary_features/test12_secondary_motif_stats_per_region_2_check.csv",
                          "test_files/motif_secondary_features/test12_secondary_motif_stats_per_file_1_check.csv",
                          "test_files/motif_secondary_features/test12_secondary_motif_stats_per_file_2_check.csv"]),
                        # count_match_list
                        ([2, 2, 2, 2])),
                        # Test 2: Multiple regions and multiple files.
                        # region_fasta_list
                        ((["test_files/motif_secondary_features/test131.fasta",
                           "test_files/motif_secondary_features/test132.fasta"]),
                        # kmer_size
                        (4),
                        # motif_list
                        (["CGCC", "CGCT", "TGCC", "TGCT"]),
                        # motif_title
                        ("YGCY"),
                        # output_file_prefix
                        ("test_files/motif_secondary_features/test13"),
                        # region_output_check_list
                        (["test_files/motif_secondary_features/test131_secondary_motif_stats_per_region_1_check.csv",
                          "test_files/motif_secondary_features/test131_secondary_motif_stats_per_region_2_check.csv",
                          "test_files/motif_secondary_features/test131_secondary_motif_stats_per_file_1_check.csv",
                          "test_files/motif_secondary_features/test131_secondary_motif_stats_per_file_2_check.csv",
                          "test_files/motif_secondary_features/test132_secondary_motif_stats_per_region_1_check.csv",
                          "test_files/motif_secondary_features/test132_secondary_motif_stats_per_region_2_check.csv",
                          "test_files/motif_secondary_features/test132_secondary_motif_stats_per_file_1_check.csv",
                          "test_files/motif_secondary_features/test132_secondary_motif_stats_per_file_2_check.csv",]),
                        # count_match_list
                        ([2, 2, 3, 3, 3, 3, 3, 3]))               
                        ])

def test_organize_files(region_fasta_list, kmer_size, motif_list, motif_title, \
                        output_file_prefix, region_output_check_list, count_match_list):
    """
        GIVEN list of motifs, a given motif size, a motif title, and a
            FASTA file of RBP binding regions that need to be analyzed.
        WHEN the function analyzes the regions, grabs all the motif
            secondary output characteristics and outputs them into
            four files.
        THEN the output files are checked for the correct data using
            information from premade files.
    """
    # Executes function and checks output.
    m_s_f.organize_files(region_fasta_list, kmer_size, motif_list, motif_title, \
                         output_file_prefix)    
    # Counts used to check output.
    line_count = 0
    # Dictionary used to save files needed to be checked.
    file_dictionary = {}
    # Grabs file name and adds them to the file_list.
    for fasta_file_name in region_fasta_list:
        # Grabs file name from region_fasta_file.
        if fasta_file_name.find("/") == -1:
            file_name = fasta_file_name
        else:    
            file_name = fasta_file_name.split("/")[-1]
        # Creates output file names.
        output_file_name_1 = f"{output_file_prefix}_{file_name}_secondary_motif_stats_per_region_1.csv"
        output_file_name_2 = f"{output_file_prefix}_{file_name}_secondary_motif_stats_per_region_2.csv"
        output_file_name_3 = (output_file_prefix + "_secondary_motif_stats_per_file_1.csv")
        output_file_name_4 = (output_file_prefix + "_secondary_motif_stats_per_file_2.csv")
        # Saves files by file name.
        file_dictionary[file_name] = [output_file_name_1, output_file_name_2, 
                                      output_file_name_3, output_file_name_4]
    # Count used to open correct check file.
    count = 0
    # Loops through and checks output of each file. Four loops because the number of
    # iterations is short.
    for file_name in file_dictionary:
        file_list = file_dictionary[file_name]
        for file in file_list:
            output_file = open(file, "r")
            # Checks each line of the output against the expected output
            # in the results_list. Counts should match.
            for line in output_file:
                # Cleans line to make it is easier to compare.
                line_clean = line.strip("\n")
                check_opened = open(region_output_check_list[count], 'r')
                for line_check in check_opened:
                    # Cleans line to make it is easier to compare.
                    line_check_clean = line_check.strip("\n")
                    if line_clean == line_check_clean:
                        line_count += 1
            # Checks line counts.
            assert count_match_list[count] == line_count
            count += 1
            line_count = 0
            output_file.close()
            check_opened.close()

"""
Tests output_stats_per_region with various nucleotide lists.

    Data Types:
    region_stats_dictionary -- Dictionary with motif and seconardary
        statistics for each region. Dictionary is organized by region.
        For nucleotide_count, nucleotide_intra_count, nucleotide_inter_count,
        see nucleotide_content for details, for motif_region_counts, see
        aggregate_motif_counts for details.
        {region: [nucleotide_count, nucleotide_intra_count,
                  nucleotide_inter_count, motif_region_counts]}.
    output_file_prefix -- Prefix for output files.
    file_name -- Name of file being analyzed.
    CSV files -- CSV file of nucleotide and motif secondary statistics.
        Contains all the statisitics in region_stats_dictionary. The titles
        of the file and an explanatation is listed below. The files are split
        into two files to make it easier to read.
        File 1- "_secondary_motif_stats_per_region_1.csv"
        Average percentages of each nucleotides surrounding motifs, within motif groups,
        and between motif groups and average length of surrounding nucleotides:
        All surrounding- All_A_%, All_C_%, All_G_%, All_T_%, 
            All_average_length
        Within motif group- Intra_A_%, Intra_C_%, Intra_G_%, Intra_T_%, 
            Intra_average_length
        Between motif group- Inter_A_%, Inter_C_%, Inter_G_%, Inter_T_%, 
            Inter_average_length
        Average nucleotide rate of each nucleotide per motif per group number:
            A_rate_per_motif_per_group, C_rate_per_motif_per_group, 
            G_rate_per_motif_per_group, T_rate_per_motif_per_group        
        Nucleotide types per group average:
            A_per_group, C_per_group, G_per_group, T_per_group
        File 2- "_secondary_motif_stats_per_region_2.csv"    
        Surrounding nucleotide type per motif type:
            A_per_CGCC, A_per_CGCT, A_per_TGCC, A_per_TGCT, C_per_CGCC, 
            C_per_CGCT, C_per_TGCC, C_per_TGCT, G_per_CGCC, G_per_CGCT, 
            G_per_TGCC, G_per_TGCT, T_per_CGCC, T_per_CGCT, T_per_TGCC, 
            T_per_TGCT
        Motif types per group average:
            CGCC_per_group, CGCT_per_group, TGCC_per_group, TGCT_per_group
        Average motif fraction per region-
            CGCC_per_region, CGCT_per_region, TGCC_per_region, TGCT_per_region    
        Average motifs per group and average surrounding nucleotides per group
            (this is the number per 100 nucleotides):
            Motif_per_group, Nucleotide_per_group
        Total average motifs per region and total average surrounding 
            nucleotides per region length:
            Total_motifs_per_region_length, Total_nucleotides_per_region_length
    region_output_check_1 -- First CSV file of expected nucleotide and motif 
        secondary statistics.
    region_output_check_2 -- Second CSV file of expected nucleotide and motif 
        secondary statistics.
    count_match_1  -- Number of lines that match between the output and the
        first check file.
    count_match_2  -- Number of lines that match between the output and the
        second check file.            
     
"""

@pytest.mark.parametrize("region_stats_dictionary, output_file_prefix, file_name, \
                          region_output_check_1, region_output_check_2, \
                          count_match_1, count_match_2", [
                        # Test 1: Single region.
                        # region_stats_dictionary
                        (({"name1": [{"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.1, "C_percentage": 20.2, 
                                      "G_percentage": 40.3, "T_percentage": 20.4,  
                                      "Total": 10, "Average_length": 10.1},
                                     {"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.0, "C_percentage": 30.0, 
                                      "G_percentage": 40.0, "T_percentage": 50.0,  
                                      "Total": 10, "Average_length": 10.2},
                                     {"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.01, "C_percentage": 20.02, 
                                      "G_percentage": 40.03, "T_percentage": 20.04,  
                                      "Total": 10, "Average_length": 10.3},
                                     {'CGCC': 1.1, 'CGCT': 1.2, 'TGCC': 1.3, 'TGCT': 1.4, 'Total_motif': 4.5,
                                      'A': 2.6, 'C': 1.7, 'G': 1.8, 'T': 6.9, 'Total_nucleotide': 10.10, 
                                      'A_rate_per_motif_per_group': 0.05, 'A_per_group': 2.0, 
                                      'C_rate_per_motif_per_group': 0.0251, 'C_per_group': 1.01, 
                                      'G_rate_per_motif_per_group': 0.025, 'G_per_group': 1.02, 
                                      'T_rate_per_motif_per_group': 0.15, 'T_per_group': 6.0, 
                                      'A_per_CGCC': 2.11, 'C_per_CGCC': 1.1, 'G_per_CGCC': 1.2, 
                                      'T_per_CGCC': 6.3, 'A_per_CGCT': 2.12, 'C_per_CGCT': 1.2, 
                                      'G_per_CGCT': 1.3, 'T_per_CGCT': 6.41, 'A_per_TGCC': 2.13, 
                                      'C_per_TGCC': 1.21, 'G_per_TGCC': 1.31, 'T_per_TGCC': 6.42, 
                                      'A_per_TGCT': 2.14, 'C_per_TGCT': 1.22, 'G_per_TGCT': 1.32, 
                                      'T_per_TGCT': 6.43, 'CGCC_per_group': 1.1, 'CGCC_per_region': 0.252, 
                                      'CGCT_per_group': 1.3, 'CGCT_per_region': 0.254, 'TGCC_per_group': 1.1, 
                                      'TGCC_per_region': 0.2521, 'TGCT_per_group': 1.03, 'TGCT_per_region': 0.254, 
                                      'motif_per_group': 4.01, 'nucleotide_per_group': 10.01, 
                                      'Total_motifs_per_region_length': 4.0, 'Total_nucleotides_per_region_length': 10.0}]}),
                        # output_file_prefix
                        ("test_files/motif_secondary_features/test1"),
                        # file_name
                        ("file1"),
                        # region_output_check_1
                        ("test_files/motif_secondary_features/test1_secondary_motif_stats_per_region_1_check.csv"),
                        # region_output_check_2
                        ("test_files/motif_secondary_features/test1_secondary_motif_stats_per_region_2_check.csv"),
                        # count_match_1
                        (2),
                        # count_match_2
                        (2)),
                        # Test 2: Multiple regions.
                        # region_stats_dictionary
                        (({"name1": [{"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.1, "C_percentage": 20.2, 
                                      "G_percentage": 40.3, "T_percentage": 20.4,  
                                      "Total": 10, "Average_length": 10.1},
                                     {"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.0, "C_percentage": 30.0, 
                                      "G_percentage": 40.0, "T_percentage": 50.0,  
                                      "Total": 10, "Average_length": 10.2},
                                     {"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.01, "C_percentage": 20.02, 
                                      "G_percentage": 40.03, "T_percentage": 20.04,  
                                      "Total": 10, "Average_length": 10.3},
                                     {'CGCC': 1.1, 'CGCT': 1.2, 'TGCC': 1.3, 'TGCT': 1.4, 'Total_motif': 4.5,
                                      'A': 2.6, 'C': 1.7, 'G': 1.8, 'T': 6.9, 'Total_nucleotide': 10.10, 
                                      'A_rate_per_motif_per_group': 0.05, 'A_per_group': 2.0, 
                                      'C_rate_per_motif_per_group': 0.0251, 'C_per_group': 1.01, 
                                      'G_rate_per_motif_per_group': 0.025, 'G_per_group': 1.02, 
                                      'T_rate_per_motif_per_group': 0.15, 'T_per_group': 6.0, 
                                      'A_per_CGCC': 2.11, 'C_per_CGCC': 1.1, 'G_per_CGCC': 1.2, 
                                      'T_per_CGCC': 6.3, 'A_per_CGCT': 2.12, 'C_per_CGCT': 1.2, 
                                      'G_per_CGCT': 1.3, 'T_per_CGCT': 6.41, 'A_per_TGCC': 2.13, 
                                      'C_per_TGCC': 1.21, 'G_per_TGCC': 1.31, 'T_per_TGCC': 6.42, 
                                      'A_per_TGCT': 2.14, 'C_per_TGCT': 1.22, 'G_per_TGCT': 1.32, 
                                      'T_per_TGCT': 6.43, 'CGCC_per_group': 1.1, 'CGCC_per_region': 0.252, 
                                      'CGCT_per_group': 1.3, 'CGCT_per_region': 0.254, 'TGCC_per_group': 1.1, 
                                      'TGCC_per_region': 0.2521, 'TGCT_per_group': 1.03, 'TGCT_per_region': 0.254, 
                                      'motif_per_group': 4.01, 'nucleotide_per_group': 10.01, 
                                      'Total_motifs_per_region_length': 4.0, 'Total_nucleotides_per_region_length': 10.0}],
                        "name2": [{"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.0, "C_percentage": 20.0, 
                                      "G_percentage": 40.0, "T_percentage": 20.0,  
                                      "Total": 10, "Average_length": 10.1},
                                     {"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.0, "C_percentage": 30.0, 
                                      "G_percentage": 40.0, "T_percentage": 50.0,  
                                      "Total": 10, "Average_length": 10.2},
                                     {"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.0, "C_percentage": 20.0, 
                                      "G_percentage": 40.0, "T_percentage": 20.0,  
                                      "Total": 10, "Average_length": 10.3},
                                     {'CGCC': 1.1, 'CGCT': 1.2, 'TGCC': 1.3, 'TGCT': 1.0, 'Total_motif': 4.5,
                                      'A': 2.0, 'C': 1.0, 'G': 1.0, 'T': 6.0, 'Total_nucleotide': 10.1, 
                                      'A_rate_per_motif_per_group': 0.05, 'A_per_group': 2.0, 
                                      'C_rate_per_motif_per_group': 0.025, 'C_per_group': 1.0, 
                                      'G_rate_per_motif_per_group': 0.025, 'G_per_group': 1.0, 
                                      'T_rate_per_motif_per_group': 0.15, 'T_per_group': 6.0, 
                                      'A_per_CGCC': 2.1, 'C_per_CGCC': 1.1, 'G_per_CGCC': 1.2, 
                                      'T_per_CGCC': 6.3, 'A_per_CGCT': 2.1, 'C_per_CGCT': 1.2, 
                                      'G_per_CGCT': 1.3, 'T_per_CGCT': 6.4, 'A_per_TGCC': 2.1, 
                                      'C_per_TGCC': 1.2, 'G_per_TGCC': 1.3, 'T_per_TGCC': 6.4, 
                                      'A_per_TGCT': 2.1, 'C_per_TGCT': 1.2, 'G_per_TGCT': 1.3, 
                                      'T_per_TGCT': 6.4, 'CGCC_per_group': 1.1, 'CGCC_per_region': 0.25, 
                                      'CGCT_per_group': 1.3, 'CGCT_per_region': 0.25, 'TGCC_per_group': 1.0, 
                                      'TGCC_per_region': 0.25, 'TGCT_per_group': 1.0, 'TGCT_per_region': 0.25, 
                                      'motif_per_group': 4.0, 'nucleotide_per_group': 10.0, 
                                      'Total_motifs_per_region_length': 4.0, 'Total_nucleotides_per_region_length': 10.0}]}),
                                  # output_file_prefix
                                  ("test_files/motif_secondary_features/test2"),
                                  # file_name
                                  ("file2"),
                                  # region_output_check_1
                                  ("test_files/motif_secondary_features/test2_secondary_motif_stats_per_region_1_check.csv"),
                                  # region_output_check_2
                                 ("test_files/motif_secondary_features/test2_secondary_motif_stats_per_region_2_check.csv"),
                                 # count_match_1
                                 (3),
                                 # count_match_2
                                 (3)),    
                        ])

def test_output_stats_per_region(region_stats_dictionary, output_file_prefix, file_name,
                                 region_output_check_1, region_output_check_2, 
                                 count_match_1, count_match_2):
    """
        GIVEN a region_stats_dictionary created by get_motif_counts, that needs
            to be written to file.
        WHEN the function extracts the information from the dictionary and
            writes it to a CSV file.
        THEN the output files are checked for the correct data and titles.
    """
    # Executes function.
    m_s_f.output_stats_per_region(region_stats_dictionary, output_file_prefix, file_name)
    # Counts used to check output.
    line_count_1 = 0
    line_count_2 = 0
    # Creates output file names.
    output_file_name_1 = f"{output_file_prefix}_{file_name}_secondary_motif_stats_per_region_1.csv"
    output_file_name_2 = f"{output_file_prefix}_{file_name}_secondary_motif_stats_per_region_2.csv"
    # Creates output files.
    output_file_1 = open(output_file_name_1, "r")
    output_file_2 = open(output_file_name_2, "r")
    # Checks each line of the output against the expected output
    # in the results_list. Counts should match.
    for line in output_file_1:
        # Cleans line to make it is easier to compare.
        line_clean = line.strip("\n")
        check1_opened = open(region_output_check_1, 'r')
        for line_check in check1_opened:
            # Cleans line to make it is easier to compare.
            line_check_clean = line_check.strip("\n")
            if line_clean == line_check_clean:
                line_count_1 += 1
    # Checks line counts.
    assert count_match_1 == line_count_1
    output_file_1.close()
    check1_opened.close()
    # Checks each line of the output against the expected output
    # in the results_list. Counts should match.
    for line in output_file_2:
        # Cleans line to make it is easier to compare.
        line_clean = line.strip("\n")
        check2_opened = open(region_output_check_2, 'r')
        for line_check in check2_opened:
            # Cleans line to make it is easier to compare.
            line_check_clean = line_check.strip("\n")
            if line_clean == line_check_clean:
                line_count_2 += 1
    # Checks line counts.
    assert count_match_2 == line_count_2
    output_file_2.close()
    check1_opened.close()   

"""
Tests output_stats_per_file with various nucleotide lists.

    Data Types:
    file_stats_dictionary_all -- Dictionary with motif and seconardary
        statistics for each input file.
    output_file_prefix -- Prefix for output files.
    CSV files -- CSV file of nucleotide and motif secondary statistics.
        Contains all the statisitics in region_stats_dictionary. The titles
        of the file and an explanatation is listed below. The files are split
        into two files to make it easier to read.
        File 1- "_secondary_motif_stats_per_region_1.csv"
        Average percentages of each nucleotides surrounding motifs, within motif groups,
        and between motif groups and average length of surrounding nucleotides:
        All surrounding- All_A_%, All_C_%, All_G_%, All_T_%, 
            All_average_length
        Within motif group- Intra_A_%, Intra_C_%, Intra_G_%, Intra_T_%, 
            Intra_average_length
        Between motif group- Inter_A_%, Inter_C_%, Inter_G_%, Inter_T_%, 
            Inter_average_length
        Average nucleotide rate of each nucleotide per motif per group number:
            A_rate_per_motif_per_group, C_rate_per_motif_per_group, 
            G_rate_per_motif_per_group, T_rate_per_motif_per_group        
        Nucleotide types per group average:
            A_per_group, C_per_group, G_per_group, T_per_group
        File 2- "_secondary_motif_stats_per_region_2.csv"    
        Surrounding nucleotide type per motif type:
            A_per_CGCC, A_per_CGCT, A_per_TGCC, A_per_TGCT, C_per_CGCC, 
            C_per_CGCT, C_per_TGCC, C_per_TGCT, G_per_CGCC, G_per_CGCT, 
            G_per_TGCC, G_per_TGCT, T_per_CGCC, T_per_CGCT, T_per_TGCC, 
            T_per_TGCT
        Motif types per group average:
            CGCC_per_group, CGCT_per_group, TGCC_per_group, TGCT_per_group
        Average motif fraction per region-
            CGCC_per_region, CGCT_per_region, TGCC_per_region, TGCT_per_region    
        Average motifs per group and average surrounding nucleotides per group
            (this is the number per 100 nucleotides):
            Motif_per_group, Nucleotide_per_group
        Total average motifs per region and total average surrounding 
            nucleotides per region length:
            Total_motifs_per_region_length, Total_nucleotides_per_region_length
        file_output_check_1 -- First CSV file of expected nucleotide and motif 
            secondary statistics.
        file_output_check_2 -- Secon CSV file of expected nucleotide and motif 
            secondary statistics.
        count_match_1  -- Number of lines that match between the output and the
            first check file.
        count_match_2  -- Number of lines that match between the output and the
            second check file.     
     
"""

@pytest.mark.parametrize("file_stats_dictionary_all, output_file_prefix, \
                          file_output_check_1, file_output_check_2, \
                          count_match_1, count_match_2", [
                          # Test 1: Single region with two motif groups. One
                          # group with two motifs, the other with one.
                          # region_stats_dictionary_all
                          (({"name1": [{"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.1, "C_percentage": 20.2, 
                                      "G_percentage": 40.3, "T_percentage": 20.4,  
                                      "Total": 10, "Average_length": 10.1},
                                     {"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.0, "C_percentage": 30.0, 
                                      "G_percentage": 40.0, "T_percentage": 50.0,  
                                      "Total": 10, "Average_length": 10.2},
                                     {"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.01, "C_percentage": 20.02, 
                                      "G_percentage": 40.03, "T_percentage": 20.04,  
                                      "Total": 10, "Average_length": 10.3},
                                     {'CGCC': 1.1, 'CGCT': 1.2, 'TGCC': 1.3, 'TGCT': 1.4, 'Total_motif': 4.5,
                                      'A': 2.6, 'C': 1.7, 'G': 1.8, 'T': 6.9, 'Total_nucleotide': 10.10, 
                                      'A_rate_per_motif_per_group': 0.05, 'A_per_group': 2.0, 
                                      'C_rate_per_motif_per_group': 0.0251, 'C_per_group': 1.01, 
                                      'G_rate_per_motif_per_group': 0.025, 'G_per_group': 1.02, 
                                      'T_rate_per_motif_per_group': 0.15, 'T_per_group': 6.0, 
                                      'A_per_CGCC': 2.11, 'C_per_CGCC': 1.1, 'G_per_CGCC': 1.2, 
                                      'T_per_CGCC': 6.3, 'A_per_CGCT': 2.12, 'C_per_CGCT': 1.2, 
                                      'G_per_CGCT': 1.3, 'T_per_CGCT': 6.41, 'A_per_TGCC': 2.13, 
                                      'C_per_TGCC': 1.21, 'G_per_TGCC': 1.31, 'T_per_TGCC': 6.42, 
                                      'A_per_TGCT': 2.14, 'C_per_TGCT': 1.22, 'G_per_TGCT': 1.32, 
                                      'T_per_TGCT': 6.43, 'CGCC_per_group': 1.1, 'CGCC_per_region': 0.252, 
                                      'CGCT_per_group': 1.3, 'CGCT_per_region': 0.254, 'TGCC_per_group': 1.1, 
                                      'TGCC_per_region': 0.2521, 'TGCT_per_group': 1.03, 'TGCT_per_region': 0.254, 
                                      'motif_per_group': 4.01, 'nucleotide_per_group': 10.01, 
                                      'Total_motifs_per_region_length': 4.0, 'Total_nucleotides_per_region_length': 10.0}]}),
                        # output_file_prefix
                        ("test_files/motif_secondary_features/test1"),
                        # region_output_check_1
                        ("test_files/motif_secondary_features/test1_secondary_motif_stats_per_file_1_check.csv"),
                        # region_output_check_2
                        ("test_files/motif_secondary_features/test1_secondary_motif_stats_per_file_2_check.csv"),
                        # count_match_1
                        (2),
                        # count_match_2
                        (2)),
                        # Test 2: Multiple regions.
                        # region_stats_dictionary_all
                        (({"name1": [{"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.1, "C_percentage": 20.2, 
                                      "G_percentage": 40.3, "T_percentage": 20.4,  
                                      "Total": 10, "Average_length": 10.1},
                                     {"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.0, "C_percentage": 30.0, 
                                      "G_percentage": 40.0, "T_percentage": 50.0,  
                                      "Total": 10, "Average_length": 10.2},
                                     {"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.01, "C_percentage": 20.02, 
                                      "G_percentage": 40.03, "T_percentage": 20.04,  
                                      "Total": 10, "Average_length": 10.3},
                                     {'CGCC': 1.1, 'CGCT': 1.2, 'TGCC': 1.3, 'TGCT': 1.4, 'Total_motif': 4.5,
                                      'A': 2.6, 'C': 1.7, 'G': 1.8, 'T': 6.9, 'Total_nucleotide': 10.10, 
                                      'A_rate_per_motif_per_group': 0.05, 'A_per_group': 2.0, 
                                      'C_rate_per_motif_per_group': 0.0251, 'C_per_group': 1.01, 
                                      'G_rate_per_motif_per_group': 0.025, 'G_per_group': 1.02, 
                                      'T_rate_per_motif_per_group': 0.15, 'T_per_group': 6.0, 
                                      'A_per_CGCC': 2.11, 'C_per_CGCC': 1.1, 'G_per_CGCC': 1.2, 
                                      'T_per_CGCC': 6.3, 'A_per_CGCT': 2.12, 'C_per_CGCT': 1.2, 
                                      'G_per_CGCT': 1.3, 'T_per_CGCT': 6.41, 'A_per_TGCC': 2.13, 
                                      'C_per_TGCC': 1.21, 'G_per_TGCC': 1.31, 'T_per_TGCC': 6.42, 
                                      'A_per_TGCT': 2.14, 'C_per_TGCT': 1.22, 'G_per_TGCT': 1.32, 
                                      'T_per_TGCT': 6.43, 'CGCC_per_group': 1.1, 'CGCC_per_region': 0.252, 
                                      'CGCT_per_group': 1.3, 'CGCT_per_region': 0.254, 'TGCC_per_group': 1.1, 
                                      'TGCC_per_region': 0.2521, 'TGCT_per_group': 1.03, 'TGCT_per_region': 0.254, 
                                      'motif_per_group': 4.01, 'nucleotide_per_group': 10.01, 
                                      'Total_motifs_per_region_length': 4.0, 'Total_nucleotides_per_region_length': 10.0}],
                           "name2": [{"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.0, "C_percentage": 20.0, 
                                      "G_percentage": 40.0, "T_percentage": 20.0,  
                                      "Total": 10, "Average_length": 10.1},
                                     {"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.0, "C_percentage": 30.0, 
                                      "G_percentage": 40.0, "T_percentage": 50.0,  
                                      "Total": 10, "Average_length": 10.2},
                                     {"A": 2, "C": 2, "G": 4, "T": 2,
                                      "A_percentage": 20.0, "C_percentage": 20.0, 
                                      "G_percentage": 40.0, "T_percentage": 20.0,  
                                      "Total": 10, "Average_length": 10.3},
                                     {'CGCC': 1.1, 'CGCT': 1.2, 'TGCC': 1.3, 'TGCT': 1.0, 'Total_motif': 4.5,
                                      'A': 2.0, 'C': 1.0, 'G': 1.0, 'T': 6.0, 'Total_nucleotide': 10.1, 
                                      'A_rate_per_motif_per_group': 0.05, 'A_per_group': 2.0, 
                                      'C_rate_per_motif_per_group': 0.025, 'C_per_group': 1.0, 
                                      'G_rate_per_motif_per_group': 0.025, 'G_per_group': 1.0, 
                                      'T_rate_per_motif_per_group': 0.15, 'T_per_group': 6.0, 
                                      'A_per_CGCC': 2.1, 'C_per_CGCC': 1.1, 'G_per_CGCC': 1.2, 
                                      'T_per_CGCC': 6.3, 'A_per_CGCT': 2.1, 'C_per_CGCT': 1.2, 
                                      'G_per_CGCT': 1.3, 'T_per_CGCT': 6.4, 'A_per_TGCC': 2.1, 
                                      'C_per_TGCC': 1.2, 'G_per_TGCC': 1.3, 'T_per_TGCC': 6.4, 
                                      'A_per_TGCT': 2.1, 'C_per_TGCT': 1.2, 'G_per_TGCT': 1.3, 
                                      'T_per_TGCT': 6.4, 'CGCC_per_group': 1.1, 'CGCC_per_region': 0.25, 
                                      'CGCT_per_group': 1.3, 'CGCT_per_region': 0.25, 'TGCC_per_group': 1.0, 
                                      'TGCC_per_region': 0.25, 'TGCT_per_group': 1.0, 'TGCT_per_region': 0.25, 
                                      'motif_per_group': 4.0, 'nucleotide_per_group': 10.0, 
                                      'Total_motifs_per_region_length': 4.0, 'Total_nucleotides_per_region_length': 10.0}]}),
                        # output_file_prefix
                        ("test_files/motif_secondary_features/test2"),
                        # region_output_check_1
                        ("test_files/motif_secondary_features/test2_secondary_motif_stats_per_file_1_check.csv"),
                        # region_output_check_2
                        ("test_files/motif_secondary_features/test2_secondary_motif_stats_per_file_2_check.csv"),
                        # count_match_1
                        (3),
                        # count_match_2
                        (3)),                                   
                        ])

def test_output_stats_per_file(file_stats_dictionary_all, output_file_prefix,
                               file_output_check_1, file_output_check_2,
                               count_match_1, count_match_2):
    """
        GIVEN list of motifs, a given motif size, a motif title, and a
            FASTA file of RBP binding regions that need to be analyzed.
        WHEN the function organizes the nucleotides surrounding each motif
            and the motifs themselves into four lists.
        THEN the output lists are checked for the correct nucleotides and
            motifs.
    """
     # Executes function.
    m_s_f.output_stats_per_file(file_stats_dictionary_all, output_file_prefix)
    # Counts used to check output.
    line_count_1 = 0
    line_count_2 = 0
    # Creates output file names.
    output_file_name_1 = (output_file_prefix + "_secondary_motif_stats_per_file_1.csv")
    output_file_name_2 = (output_file_prefix + "_secondary_motif_stats_per_file_2.csv")
    # Creates output files.
    output_file_1 = open(output_file_name_1, "r")
    output_file_2 = open(output_file_name_2, "r")
    # Checks each line of the output against the expected output
    # in the results_list. Counts should match.
    for line in output_file_1:
        # Cleans line to make it is easier to compare.
        line_clean = line.strip("\n")
        check1_opened = open(file_output_check_1, 'r')
        for line_check in check1_opened:
            # Cleans line to make it is easier to compare.
            line_check_clean = line_check.strip("\n")
            if line_clean == line_check_clean:
                line_count_1 += 1
    # Checks line counts.
    assert count_match_1 == line_count_1
    output_file_1.close()
    check1_opened.close()
    # Checks each line of the output against the expected output
    # in the results_list. Counts should match.
    for line in output_file_2:
        # Cleans line to make it is easier to compare.
        line_clean = line.strip("\n")
        check2_opened = open(file_output_check_2, 'r')
        for line_check in check2_opened:
            # Cleans line to make it is easier to compare.
            line_check_clean = line_check.strip("\n")
            if line_clean == line_check_clean:
                line_count_2 += 1
    # Checks line counts.
    assert count_match_2 == line_count_2
    output_file_2.close()
    check1_opened.close()            