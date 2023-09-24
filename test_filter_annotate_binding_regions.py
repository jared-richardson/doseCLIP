"""test_filter_annotate_binding_regions.py- Unit tests for
    filter_annotate_binding_regions.py
"""
import pytest
# Imports all functions from join_binding_regions.py.
import filter_annotate_binding_regions as f_a_b_r

"""
Tests add_events_to_dictionary() with different DeSeq2 annotation files and GTF files.

Data Types:
    filtered_file -- DeSeq2 produced SM counts file or CLIP file (CSV) for a single protein
        concentration (produced with replicates). This should have been analyzed
        and normalized with DeSeq2 with the regular CLIP samples and the the SM 
        Input samples. Also, this file could be the filtered differential expression
        output file or the normalized counts file.
    regular_gtf_file -- Gencode/Ensembl/UCSC provided GTF file. Should be the GTF
        annotation for the target alignment sequence. The file is used to annotate
        the binding regions.
    region_dictionary -- Filled output dictionary with each binding region names
        as the key and the clean line as the value.
        {"region_name": "line_clean=gene_line=sub_gene_line"}.
    title_line -- String containing the title line of file. Used for output.      
    """

@pytest.mark.parametrize("filtered_file, regular_gtf_file, region_dictionary, \
                         title_line", [
                         # Test 1: One matching event with a normalized counts file
                         # in an exon region.
                         # filtered_file
                         (('test_files/filter_annotate_binding_regions/sm_test1_normalized_counts.csv'),
                           # regular_gtf_file
                           ('test_files/filter_annotate_binding_regions/gtf_test1.csv'),
                          # region_dictionary
                          ({'chr1~100~300~+': 'chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik=exon'}),
                          # title_line
                          ',x50_sample_1,x50_sample_2,x50_sample_3,x5_sample_1,x5_sample_2,x5_sample_3'),
                         # Test 2: One matching event with a normalized counts file
                         # in an exon region. Contains multiple additional lines in
                         # DeSeq2 file and counts file that do not match.
                         # filtered_file
                         (('test_files/filter_annotate_binding_regions/sm_test2_normalized_counts.csv'),
                           # regular_gtf_file
                           ('test_files/filter_annotate_binding_regions/gtf_test2.csv'),
                          # region_dictionary
                          ({'chr1~100~300~+': 'chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik=exon',
                            'chr1~100~300~-': 'chr1~100~300~-,100.0,100.3,100.2,50.1,50.2,50.4==',
                            'chr1~3100~3300~+': 'chr1~3100~3300~+,100.0,100.3,100.2,50.1,50.2,50.4==',
                            'chr2~100~300~+': 'chr2~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4=='}),
                          # title_line
                          ',x50_sample_1,x50_sample_2,x50_sample_3,x5_sample_1,x5_sample_2,x5_sample_3'),
                         # Test 3: Multiple matching event with a normalized counts file
                         # in multiple, single annotated regions. Contains other additional lines in
                         # DeSeq2 file and counts file that do not match. Contains 5' UTR.
                         # filtered_file
                         (('test_files/filter_annotate_binding_regions/sm_test2_normalized_counts.csv'),
                           # regular_gtf_file
                           ('test_files/filter_annotate_binding_regions/gtf_test3.csv'),
                          # region_dictionary
                          ({'chr1~100~300~+': 'chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik=exon',
                            'chr1~100~300~-': 'chr1~100~300~-,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000089699.2~lncRNA~Gm1992=',
                            'chr1~3100~3300~+': 'chr1~3100~3300~+,100.0,100.3,100.2,50.1,50.2,50.4==',
                            'chr2~100~300~+': 'chr2~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693~snRNA~gene=5_UTR'}),
                          # title_line
                          ',x50_sample_1,x50_sample_2,x50_sample_3,x5_sample_1,x5_sample_2,x5_sample_3'),
                         # Test 4: Multiple matching event with a normalized counts file
                         # in multiple, multiple overlapping annotated regions. Some annotated regions
                         # are smaller and/or larger than binding region. Contains multiple 
                         # additional lines in DeSeq2 file and counts file that do not match. Contains
                         # 5' UTR and 3' UTR.
                         # filtered_file
                         (('test_files/filter_annotate_binding_regions/sm_test3_normalized_counts.csv'),
                           # regular_gtf_file
                           ('test_files/filter_annotate_binding_regions/gtf_test4.csv'),
                          # region_dictionary
                          ({'chr1~50~300~+': 'chr1~50~300~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=5_UTR~exon~3_UTR',
                            'chr1~100~300~-': 'chr1~100~300~-,100.0,100.3,100.2,50.1,50.2,50.4==',
                            'chr1~3100~3300~+': 'chr1~3100~3300~+,100.0,100.3,100.2,50.1,50.2,50.4==',
                            'chr2~4599~4640~+': 'chr2~4599~4640~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693~snRNA~gene=3_UTR'}),
                          # title_line
                          ',x50_sample_1,x50_sample_2,x50_sample_3,x5_sample_1,x5_sample_2,x5_sample_3'),
                         # Test 5: Same test as four but with DeSeq2 differential expression file, instead
                         # of normalized counts. Multiple matching event with a normalized counts file
                         # in multiple, multiple overlapping annotated regions. Some annotated regions
                         # are smaller and/or larger than binding region. Contains multiple 
                         # additional lines in DeSeq2 file and counts file that do not match. Contains
                         # 5' UTR and 3' UTR.
                         # filtered_file
                         (('test_files/filter_annotate_binding_regions/sm_test3_differential_expression.csv'),
                           # regular_gtf_file
                           ('test_files/filter_annotate_binding_regions/gtf_test4.csv'),
                          # region_dictionary
                          ({'chr1~50~300~+': 'chr1~50~300~+,200,5,1,4,0.005,0.05=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=5_UTR~exon~3_UTR',
                            'chr1~100~300~-': 'chr1~100~300~-,200,5,1,4,0.005,0.05==',
                            'chr1~3100~3300~+': 'chr1~3100~3300~+,200,5,1,4,0.005,0.05==',
                            'chr2~4599~4640~+': 'chr2~4599~4640~+,200,5,1,4,0.005,0.05=ENSMUSG00000102693~snRNA~gene=3_UTR'}),
                          # title_line
                          ",baseMean,log2FoldChange,lfcSE,stat,pvalue,padj")
                         ])

def test_add_events_to_dictionary(filtered_file, regular_gtf_file,
                                  region_dictionary, title_line):
    """
        GIVEN a DeSeq2 file with binding regions that needs to be annotated.
        WHEN the function finds overlapping genes and sub-gene types that
            overlap with the binding region and matches. The information
            is then added to a dictionary.
        THEN the output dictionary is checked for the correct region and 
            gene information. The title line string from the input file 
            is also checked for correctness.
    """
    assert f_a_b_r.add_events_to_dictionary(filtered_file, regular_gtf_file) == (region_dictionary, title_line)

"""
Tests output_file() with different combinations of DeSeq2 input files and GTF annotations.

Data Types:
    sm_filtered_file -- DeSeq2 produced SM counts file (CSV) for a single protein
        concentration (produced with replicates). This should have been analyzed
        and normalized with DeSeq2 with the regular CLIP samples and the the SM 
        Input samples. Also, this file could be the filtered differential expression
        output file or the normalized counts file. 
    clip_regular_file -- DeSeq2 produced CLIP counts file (CSV) for a single protein
        concentration (produced with replicates). This file should have been analyzed
        and normalized with DeSeq2 without the SM Input samples. Also, this file could 
        be the filtered differential expression output file or the normalized
        counts file. These are the numerical values that are primarily used for later
        analyses. This is an OPTIONAL parameter, as the script is able to just annotate 
        the SM filtered counts as well.
    regular_gtf_file -- Gencode/Ensembl/UCSC provided GTF file. Should be the GTF
        annotation for the target alignment sequence. The file is used to annotate
        the binding regions.
    sample_keyword -- Keyword to tell the function to process the sm_filtered_file or
        the clip_regular_file. Keyword options are "SM" for sm_filtered_file and
        "CLIP" for clip_regular_file.
    sm_filtered_file_out -- Output normalized SM filtered counts or SM differential 
        expression file (this is output from sm_filtered_file). The file will 
        be annotated. The file is in CSV format.
    clip_regular_file_out -- Output normalized regular CLIP filtered counts or 
        CLIP differential expression file (this is output from clip_regular_file). 
        The file will be annotated and filtered with the SM filtered binding regions.
        The file is in CSV format. This is an OPTIONAL parameter, as the script is 
        able to just annotate the SM filtered counts as well.
    sm_filtered_file_out_count -- Sub-gene count output normalized SM filtered counts 
        or SM differential expression file (this is output from sm_filtered_file). 
        The file will be annotated. The file is in CSV format.
    clip_regular_file_out_count -- Sub-gene count output normalized regular CLIP 
        filtered counts or CLIP differential expression file (this is output from 
        clip_regular_file). The file will be annotated and filtered with the SM 
        filtered binding regions. The file is in CSV format. This is an OPTIONAL 
        parameter, as the script is able to just annotate the SM filtered counts 
        as well.    
    results_list -- List of results expected in function produced files for the
        SM filtered samples.
    count_match -- Expected count of matching lines in output and premade
        testing files for the SM filtered samples.
    results_list_count -- List of results expected in function for produced sub-gene
        counts file for the SM filtered samples.
    count_match_count -- Expected count of matching lines in output and premade
        testing files for produced sub-gene counts file for the SM filtered samples.
    results_list_clip -- List of results expected in function produced files for the
        CLIP output.
    count_match_clip -- Expected count of matching lines in output and premade
        testing files for the CLIP output.
    results_list_clip_count -- List of results expected in function produced files for the
        CLIP sub-gene counts output.
    count_match_clip_count -- Expected count of matching lines in output and premade
        testing files for the CLIP sub-gene counts output.
    """

@pytest.mark.parametrize("sm_filtered_file, regular_gtf_file, sm_filtered_file_out, sample_keyword, \
                          clip_regular_file, clip_regular_file_out, sm_filtered_file_out_count, \
                          clip_regular_file_out_count, results_list, count_match, results_list_count, \
                          count_match_count, results_list_clip, count_match_clip, results_list_clip_count, \
                          count_match_clip_count", [
                         # Test 1: One matching event with a normalized counts file, no regular CLIP files.
                         # in an exon region.
                         # sm_filtered_file
                         (("test_files/filter_annotate_binding_regions/sm_test1_normalized_counts.csv"),
                           # regular_gtf_file
                           ("test_files/filter_annotate_binding_regions/gtf_test1.csv"),
                           # sm_filtered_file_out
                           ("test_files/filter_annotate_binding_regions/sm_test1_normalized_counts_out.csv"),
                           # sample_keyword
                           "SM",
                           # clip_regular_file
                           None,
                           # clip_regular_file_out
                           None,
                           # sm_filtered_file_out_count
                           ("test_files/filter_annotate_binding_regions/sm_test1_normalized_counts_out.annotation_counts.csv"),
                           # clip_regular_file_out_count
                           (""),
                           # results_list
                           [",x50_sample_1,x50_sample_2,x50_sample_3,x5_sample_1,x5_sample_2,x5_sample_3,gene_id,gene_type,gene_name,sub_gene_type,all_genes,all_sub_gene_types",
                            "chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693.2,TEC,4933401J01Rik,exon,,"],
                           # count_match
                           (2),
                           # results_list_count
                           ["Sub-Gene Type,Sub-Gene Count", "exon,1"],
                           # count_match_count
                           (2), 
                           # results_list_clip
                           [],
                           # count_match_clip
                           (0),
                           # results_list_clip_count
                           [],
                           # count_match_clip_count
                           (0)),
                         # Test 2: Multiple matching events with multiple genes/sub-genes with a 
                         # normalized counts file, no regular CLIP files.
                         # in an exon region.
                         # sm_filtered_file
                         (("test_files/filter_annotate_binding_regions/sm_test2_normalized_counts.csv"),
                           # regular_gtf_file
                           ("test_files/filter_annotate_binding_regions/gtf_test4.csv"),
                           # sm_filtered_file_out
                           ("test_files/filter_annotate_binding_regions/sm_test2_normalized_counts_out.csv"),
                           # sample_keyword
                           "SM",
                           # clip_regular_file
                           None,
                           # clip_regular_file_out
                           None,
                           # sm_filtered_file_out_count
                           ("test_files/filter_annotate_binding_regions/sm_test2_normalized_counts_out.annotation_counts.csv"),
                           # clip_regular_file_out_count
                           (""),
                           # results_list
                           [",x50_sample_1,x50_sample_2,x50_sample_3,x5_sample_1,x5_sample_2,x5_sample_3,gene_id,gene_type,gene_name,sub_gene_type,all_genes,all_sub_gene_types",
                            "chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693.2,TEC,4933401J01Rik,exon,ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992,exon~3_UTR",
                            "chr1~100~300~-,100.0,100.3,100.2,50.1,50.2,50.4,,,,,,",
                            "chr1~3100~3300~+,100.0,100.3,100.2,50.1,50.2,50.4,,,,,,",
                            "chr2~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693,snRNA,gene,,,"],
                           # count_match
                           (5),
                           # results_list_count
                           ["Sub-Gene Type,Sub-Gene Count", "3_UTR,1", "intron,1", "intergenic,2"],
                           # count_match_count
                           (4),
                           # results_list_clip
                           [],
                           # count_match_clip
                           (0),
                           # results_list_clip_count
                           [],
                           # count_match_clip_count
                           (0)),
                         # Test 3: One matching event with a normalized counts file, with regular CLIP files.
                         # in an exon region.
                         # sm_filtered_file
                         (("test_files/filter_annotate_binding_regions/sm_test1_normalized_counts.csv"),
                           # regular_gtf_file
                           ("test_files/filter_annotate_binding_regions/gtf_test1.csv"),
                           # sm_filtered_file_out
                           ("test_files/filter_annotate_binding_regions/sm_test1_normalized_counts_out.csv"),
                           # sample_keyword
                           "CLIP",
                           # clip_regular_file
                           "test_files/filter_annotate_binding_regions/clip_test1_normalized_counts.csv",
                           # clip_regular_file_out
                           "test_files/filter_annotate_binding_regions/clip_test1_normalized_counts_out.csv",
                           # sm_filtered_file_out_count
                           ("test_files/filter_annotate_binding_regions/sm_test1_normalized_counts_out.annotation_counts.csv"),
                           # clip_regular_file_out_count
                           ("test_files/filter_annotate_binding_regions/clip_test1_normalized_counts_out.annotation_counts.csv"),
                           # results_list
                           [",x50_sample_1,x50_sample_2,x50_sample_3,x5_sample_1,x5_sample_2,x5_sample_3,gene_id,gene_type,gene_name,sub_gene_type,all_genes,all_sub_gene_types",
                            "chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693.2,TEC,4933401J01Rik,exon,,"],
                           # count_match
                           (2),
                           # results_list_count
                           ["Sub-Gene Type,Sub-Gene Count", "exon,1"],
                           # count_match_count
                           (2),
                           # results_list_clip
                           [",x50_sample_1,x50_sample_2,x50_sample_3,x5_sample_1,x5_sample_2,x5_sample_3,gene_id,gene_type,gene_name,sub_gene_type,all_genes,all_sub_gene_types",
                            "chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693.2,TEC,4933401J01Rik,exon,,"],
                           # count_match_clip
                           (2),
                           # results_list_clip_count
                           ["Sub-Gene Type,Sub-Gene Count", "exon,1"],
                           # count_match_clip_count
                           (2)),
                         # Test 4: Multiple matching events with multiple genes/sub-genes with a 
                         # normalized counts file, with regular CLIP files.
                         # in an exon region.
                         # sm_filtered_file
                         (("test_files/filter_annotate_binding_regions/sm_test4_normalized_counts.csv"),
                           # regular_gtf_file
                           ("test_files/filter_annotate_binding_regions/gtf_test4.csv"),
                           # sm_filtered_file_out
                           ("test_files/filter_annotate_binding_regions/sm_test2_normalized_counts_out.csv"),
                           # sample_keyword
                           "CLIP",
                           # clip_regular_file
                           "test_files/filter_annotate_binding_regions/clip_test2_normalized_counts.csv",
                           # clip_regular_file_out
                           "test_files/filter_annotate_binding_regions/clip_test2_normalized_counts_out.csv",
                           # sm_filtered_file_out_count
                           ("test_files/filter_annotate_binding_regions/sm_test2_normalized_counts_out.annotation_counts.csv"),
                           # clip_regular_file_out_count
                           ("test_files/filter_annotate_binding_regions/clip_test2_normalized_counts_out.annotation_counts.csv"),
                           # results_list
                           [",x50_sample_1,x50_sample_2,x50_sample_3,x5_sample_1,x5_sample_2,x5_sample_3,gene_id,gene_type,gene_name,sub_gene_type,all_genes,all_sub_gene_types",
                            "chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693.2,TEC,4933401J01Rik,exon,ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992,exon~3_UTR",
                            "chr1~100~300~-,100.0,100.3,100.2,50.1,50.2,50.4,,,,,,",
                            "chr1~3100~3300~+,100.0,100.3,100.2,50.1,50.2,50.4,,,,,,",
                            "chr2~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693,snRNA,gene,,,"],
                           # count_match
                           (5),
                           # results_list_count
                           ["Sub-Gene Type,Sub-Gene Count", "3_UTR,1", "intron,1", "intergenic,2"],
                           # count_match_count
                           (4),
                           # results_list_clip
                           [",x50_sample_1,x50_sample_2,x50_sample_3,x5_sample_1,x5_sample_2,x5_sample_3,gene_id,gene_type,gene_name,sub_gene_type,all_genes,all_sub_gene_types",
                            "chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693.2,TEC,4933401J01Rik,exon,ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992,exon~3_UTR",
                            "chr1~100~300~-,100.0,100.3,100.2,50.1,50.2,50.4,,,,,,",
                            "chr1~3100~3300~+,100.0,100.3,100.2,50.1,50.2,50.4,,,,,,",
                            "chr2~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693,snRNA,gene,,,"],
                           # count_match_clip
                           (5),
                           # results_list_clip_count
                           ["Sub-Gene Type,Sub-Gene Count", "3_UTR,1", "intron,1", "intergenic,2"],
                           # count_match_clip_count
                           (4)),
                         ])

def test_output_file(sm_filtered_file, regular_gtf_file, sm_filtered_file_out, sample_keyword,
                          clip_regular_file, clip_regular_file_out, sm_filtered_file_out_count,
                          clip_regular_file_out_count, results_list, count_match, results_list_count,
                          count_match_count, results_list_clip, count_match_clip, results_list_clip_count,
                          count_match_clip_count):
    """
        GIVEN a DeSeq2 produced file that was SM filtered, with or without additional
            regular, non-filtered CLIP DeSeq2 produced file, and a GTF annotation, the
            files are annotated and the regular CLIP file is filtered using the SM filtered
            file.
        WHEN the function annotates the DeSeq2 input binding regions with the GTF features,
            and filters the addition non-filtered CLIP DeSeq2 events (if input) and outputs
            them in a file.
        THEN the output file is checked for the correct output by
            counting the number of perfectly matching lines.
    """
    # Sets line_count to 0. This is the count of how many
    # lines match in the output and the expected output file.
    line_count = 0
    # Executes function.
    f_a_b_r.output_file(sm_filtered_file, regular_gtf_file, sm_filtered_file_out, sample_keyword,
                        clip_regular_file, clip_regular_file_out)
    # Opens testing annotated SM filtered file to iterate through 
    # to check against the expected output in the results_list.
    read1_opened = open(sm_filtered_file_out, 'r')
    # Checks each line of the output against the expected output
    # in the results_list. Counts should match.
    for line in read1_opened:
        # Cleans line to make it is easier to compare.
        line_clean = line.strip("\n")
        if line_clean in results_list:
            line_count += 1
    # Checks line counts.
    assert count_match == line_count
    read1_opened.close()
    # Resets line_count.
    line_count = 0
    # Opens testing SM filtered sub-gene counts file to 
    # iterate through to check against the expected output in 
    # the results_list.
    read1_opened = open(sm_filtered_file_out_count, 'r')
    # Checks each line of the output against the expected output
    # in the results_list. Counts should match.
    for line in read1_opened:
        # Cleans line to make it is easier to compare.
        line_clean = line.strip("\n")
        if line_clean in results_list_count:
            line_count += 1
    # Checks line counts.
    assert count_match_count == line_count
    read1_opened.close()
    # Resets line_count.
    line_count = 0
    # If CLIP keyword is provided, checks file for correct output.
    if sample_keyword == "CLIP":
        read1_opened = open(clip_regular_file_out, 'r')
        # Checks each line of the output against the expected output
        # in the results_list. Counts should match.
        for line in read1_opened:
            # Cleans line to make it is easier to compare.
            line_clean = line.strip("\n")
            if line_clean in results_list_clip:
                line_count += 1
        # Checks line counts.
        assert count_match_clip == line_count
        read1_opened.close() 
        # Resets line_count.
        line_count = 0
        # Opens testing CLIP sub-gene counts file to 
        # iterate through to check against the expected output in 
        # the results_list.
        read1_opened = open(clip_regular_file_out_count, 'r')
        # Checks each line of the output against the expected output
        # in the results_list. Counts should match.
        for line in read1_opened:
            # Cleans line to make it is easier to compare.
            line_clean = line.strip("\n")
            if line_clean in results_list_clip_count:
                line_count += 1
        # Checks line counts.
        assert count_match_clip_count == line_count
        read1_opened.close()

"""
Tests count_and_output_sub_genes() with different combinations of DeSeq2 input files and GTF annotations.

Data Types:
    region_dictionary -- Filled output dictionary with each binding region names
        as the key and the clean line as the value.
        {"region_name": "line_clean=gene_line=sub_gene_line"}.
    sm_filtered_file_out -- Output of annotation counts for normalized SM filtered 
        counts or SM differential expression file (this is output from sm_filtered_file). 
        The file is in CSV format. File contains a ".annotation_counts.csv" suffix. If no
        input (both files will not be input at once), "NA" is input and output will be
        ignored.
        File Format: Sub-Gene Type, Sub-Gene Count
    clip_regular_file_out -- Output of annotation counts for normalized regular CLIP 
        filtered counts or CLIP differential expression file (this is output from 
        clip_regular_file). The file will be annotated and filtered with the SM 
        filtered binding regions. The file is in CSV format. This is an OPTIONAL 
        parameter, as the script is able to just annotate the SM filtered counts 
        as well. File contains a ".annotation_counts.csv" suffix. If no input 
        (both files will not be input at once), "NA" is input and output will be
        ignored.
        File Format: Sub-Gene Type, Sub-Gene Count
    results_list -- List of results expected in function produced files.
    count_match -- Expected count of matching lines in output and premade
        testing files                 
    """

@pytest.mark.parametrize("region_dictionary, sm_filtered_file_out, clip_regular_file_out, \
                          results_list, count_match", [
                          # Test 1: One event with a normalized counts file.
                          # region_dictionary
                          (({'chr1~100~300~+': 'chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik=exon'}),
                           # sm_filtered_file_out
                           ("test_files/filter_annotate_binding_regions/sm_test1_normalized_counts_out.csv"),
                           # clip_regular_file_out
                           None,
                           # results_list
                           ["Sub-Gene Type,Sub-Gene Count", "exon,1"],
                           # count_match
                           (2)),
                          # Test 2: Multiple events in different orders with a normalized counts file.
                          # region_dictionary
                          (({'chr1~100~300~+': 'chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik=exon',
                             'chr1~500~700~+': 'chr1~500~700~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=CDS~exon',
                             'chr1~800~900~+': 'chr1~800~900~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=exon~CDS',
                             'chr1~950~975~+': 'chr1~950~975~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=CDS',
                             'chr1~1000~1200~+': 'chr1~1000~1200~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=3_UTR',
                             'chr1~1300~1500~+': 'chr1~1300~1500~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=5_UTR',
                             'chr1~2000~2200~+': 'chr1~2000~2200~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=start_codon~exon',
                             'chr1~2300~2400~+': 'chr1~2300~2400~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=start_codon',
                             'chr1~3000~3200~+': 'chr1~3000~3200~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=exon~stop_codon',
                             'chr1~3300~3500~+': 'chr1~3300~3500~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=stop_codon',
                             'chr1~4100~4300~+': 'chr1~4100~4300~+,200,5,1,4,0.005,0.05==',
                             'chr2~4100~4300~+': 'chr2~4100~4300~+,200,5,1,4,0.005,0.05=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992='}),
                           # sm_filtered_file_out
                           ("test_files/filter_annotate_binding_regions/sm_test1_normalized_counts_out.csv"),
                           # clip_regular_file_out
                           None,
                           # results_list
                           ["Sub-Gene Type,Sub-Gene Count", "exon,3", "CDS,3", "3_UTR,1", "5_UTR,1", "intron,1", "intergenic,1",
                            "stop_codon,1","start_codon,1"],
                           # count_match
                           (9)),
                          # Test 3: Multiple events with more counts in different orders with a normalized counts file.
                          # region_dictionary
                          (({'chr1~100~300~+': 'chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik=exon',
                             'chr1~500~700~+': 'chr1~500~700~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=CDS~exon',
                             'chr1~800~900~+': 'chr1~800~900~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=exon~CDS',
                             'chr1~950~975~+': 'chr1~950~975~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=CDS',
                             'chr1~1000~1200~+': 'chr1~1000~1200~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=3_UTR',
                             'chr1~1001~1201~+': 'chr1~1001~1201~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=3_UTR',
                             'chr1~1300~1500~+': 'chr1~1300~1500~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=5_UTR',
                             'chr1~1301~1501~+': 'chr1~1301~1501~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=5_UTR',
                             'chr1~1302~1502~+': 'chr1~1302~1502~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=5_UTR',
                             'chr1~2000~2200~+': 'chr1~2000~2200~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=start_codon~exon',
                             'chr1~2000~2201~+': 'chr1~2000~2201~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=start_codon~exon',
                             'chr1~2300~2400~+': 'chr1~2300~2400~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=start_codon',
                             'chr1~2300~2401~+': 'chr1~2300~2401~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=start_codon',
                             'chr1~2300~2402~+': 'chr1~2300~2402~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=start_codon',
                             'chr1~3000~3200~+': 'chr1~3000~3200~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=exon~stop_codon',
                             'chr1~3000~3201~+': 'chr1~3000~3201~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=exon~stop_codon',
                             'chr1~3000~3202~+': 'chr1~3000~3202~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=exon~stop_codon',
                             'chr1~3300~3500~+': 'chr1~3300~3500~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=stop_codon',
                             'chr1~4100~4300~+': 'chr1~4100~4300~+,200,5,1,4,0.005,0.05==',
                             'chr1~4100~4301~+': 'chr1~4100~4301~+,200,5,1,4,0.005,0.05==',
                             'chr1~4100~4302~+': 'chr1~4100~4302~+,200,5,1,4,0.005,0.05==',
                             'chr2~4100~4300~+': 'chr2~4100~4300~+,200,5,1,4,0.005,0.05=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=',
                             'chr2~4100~4301~+': 'chr2~4100~4301~+,200,5,1,4,0.005,0.05=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=',
                             'chr2~4100~4302~+': 'chr2~4100~4302~+,200,5,1,4,0.005,0.05=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992='}),
                           # sm_filtered_file_out
                           ("test_files/filter_annotate_binding_regions/sm_test1_normalized_counts_out.csv"),
                           # clip_regular_file_out
                           None,
                           # results_list
                           ["Sub-Gene Type,Sub-Gene Count", "exon,6", "CDS,3", "3_UTR,2", "5_UTR,3", "intron,3", "intergenic,3",
                            "stop_codon,1","start_codon,3"],
                           # count_match
                           (9)),
                          # Test 4: One event with a DeSeq2 differential expression file.
                          # region_dictionary
                          (({'chr1~50~300~+': 'chr1~50~300~+,200,5,1,4,0.005,0.05=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=5_UTR~exon',
                             'chr1~51~301~+': 'chr1~51~301~+,200,5,1,4,0.005,0.05=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=5_UTR~exon',
                             'chr1~100~300~-': 'chr1~100~300~-,200,5,1,4,0.005,0.05==',
                             'chr1~3100~3300~+': 'chr1~3100~3300~+,200,5,1,4,0.005,0.05==',
                             'chr2~4599~4640~+': 'chr2~4599~4640~+,200,5,1,4,0.005,0.05=ENSMUSG00000102693~snRNA~gene=3_UTR'}),
                           # sm_filtered_file_out
                           ("test_files/filter_annotate_binding_regions/sm_test1_normalized_counts_out.csv"),
                           # clip_regular_file_out
                           None,
                           # results_list
                           ["Sub-Gene Type,Sub-Gene Count", "exon,2", "intergenic,2","3_UTR,1"],
                           # count_match
                           (4)),
                         ])

def test_count_and_output_sub_genes(region_dictionary, sm_filtered_file_out, clip_regular_file_out,
                                    results_list, count_match):
    """
        GIVEN a dictionary of annotated binding locations from a DeSeq2 produced file of
            binding regions produced from the add_events_to_dictionary function.
        WHEN the function counts the number of sub-gene regions in each dictionary
            and outputs in a counts file.
        THEN the output file is checked for the correct output by
            counting the number of perfectly matching lines.
    """
    # Sets line_count to 0. This is the count of how many
    # lines match in the output and the expected output file.
    line_count = 0
    # Executes function.
    f_a_b_r.count_and_output_sub_genes(region_dictionary, sm_filtered_file_out, clip_regular_file_out)
    # Opens testing output files to iterate through to check against
    # the expected output in the the results_list. Checks for the
    # correct output file, as two DeSeq2 produced files can be input.
    if sm_filtered_file_out != None:
        sm_filtered_file_out = sm_filtered_file_out.replace(".csv", ".annotation_counts.csv")
        read1_opened = open(sm_filtered_file_out, 'r')
    else:
        clip_regular_file_out = clip_regular_file_out.replace(".csv", ".annotation_counts.csv")
        read1_opened = open(clip_regular_file_out, 'r')   
    # Checks each line of the output against the expected output
    # in the results_list. Counts should match.
    for line in read1_opened:
        # Cleans line to make it is easier to compare.
        line_clean = line.strip("\n")
        if line_clean in results_list:
            line_count += 1
    # Checks line counts.
    assert count_match == line_count
    read1_opened.close()

"""
Tests filter_annotate_binding_regions() with different combinations of DeSeq2 input files 
    and GTF annotations.

Data Types:
    sm_filtered_file -- DeSeq2 produced SM counts file (CSV) for a single protein
        concentration (produced with replicates). This should have been analyzed
        and normalized with DeSeq2 with the regular CLIP samples and the the SM 
        Input samples. Also, this file could be the filtered differential expression
        output file or the normalized counts file. 
    clip_regular_file -- DeSeq2 produced CLIP counts file (CSV) for a single protein
        concentration (produced with replicates). This file should have been analyzed
        and normalized with DeSeq2 without the SM Input samples. Also, this file could 
        be the filtered differential expression output file or the normalized
        counts file. These are the numerical values that are primarily used for later
        analyses. This is an OPTIONAL parameter, as the script is able to just annotate 
        the SM filtered counts as well.
    regular_gtf_file -- Gencode/Ensembl/UCSC provided GTF file. Should be the GTF
        annotation for the target alignment sequence. The file is used to annotate
        the binding regions.
    sm_filtered_file_out -- Output normalized SM filtered counts or SM differential 
        expression file (this is output from sm_filtered_file). The file will 
        be annotated. The file is in CSV format.
    clip_regular_file_out -- Output normalized regular CLIP filtered counts or 
        CLIP differential expression file (this is output from clip_regular_file). 
        The file will be annotated and filtered with the SM filtered binding regions.
        The file is in CSV format. This is an OPTIONAL parameter, as the script is 
        able to just annotate the SM filtered counts as well.
    results_list -- List of results expected in function produced files.
    count_match -- Expected count of matching lines in output and premade
        testing files                 
    """

@pytest.mark.parametrize("sm_filtered_file, regular_gtf_file, sm_filtered_file_out, \
                          clip_regular_file, clip_regular_file_out, results_list, count_match", [
                         # Test 1: One matching event with a normalized counts file, no regular CLIP files.
                         # in an exon region.
                         # sm_filtered_file
                         (("test_files/filter_annotate_binding_regions/sm_test1_normalized_counts.csv"),
                           # regular_gtf_file
                           ("test_files/filter_annotate_binding_regions/gtf_test1.csv"),
                           # sm_filtered_file_out
                           ("test_files/filter_annotate_binding_regions/sm_test1_normalized_counts_out.csv"),
                           # clip_regular_file
                           None,
                           # clip_regular_file_out
                           None,
                           # results_list
                           [",x50_sample_1,x50_sample_2,x50_sample_3,x5_sample_1,x5_sample_2,x5_sample_3,gene_id,gene_type,gene_name,sub_gene_type,all_genes,all_sub_gene_types",
                            "chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693.2,TEC,4933401J01Rik,exon,,"],
                           # count_match
                           (2)),
                         # Test 2: Multiple matching events with multiple genes/sub-genes with a 
                         # normalized counts file, no regular CLIP files.
                         # in an exon region.
                         # sm_filtered_file
                         (("test_files/filter_annotate_binding_regions/sm_test2_normalized_counts.csv"),
                           # regular_gtf_file
                           ("test_files/filter_annotate_binding_regions/gtf_test4.csv"),
                           # sm_filtered_file_out
                           ("test_files/filter_annotate_binding_regions/sm_test2_normalized_counts_out.csv"),
                           # clip_regular_file
                           None,
                           # clip_regular_file_out
                           None,
                           # results_list
                           [",x50_sample_1,x50_sample_2,x50_sample_3,x5_sample_1,x5_sample_2,x5_sample_3,gene_id,gene_type,gene_name,sub_gene_type,all_genes,all_sub_gene_types",
                            "chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693.2,TEC,4933401J01Rik,exon,ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992,exon~3_UTR",
                            "chr1~100~300~-,100.0,100.3,100.2,50.1,50.2,50.4,,,,,,",
                            "chr1~3100~3300~+,100.0,100.3,100.2,50.1,50.2,50.4,,,,,,",
                            "chr2~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693,snRNA,gene,,,"],
                           # count_match
                           (5)),
                         # Test 3: One matching event with a normalized counts file, with regular CLIP files.
                         # in an exon region.
                         # sm_filtered_file
                         (("test_files/filter_annotate_binding_regions/sm_test1_normalized_counts.csv"),
                           # regular_gtf_file
                           ("test_files/filter_annotate_binding_regions/gtf_test1.csv"),
                           # sm_filtered_file_out
                           ("test_files/filter_annotate_binding_regions/sm_test1_normalized_counts_out.csv"),
                           # clip_regular_file
                           "test_files/filter_annotate_binding_regions/clip_test1_normalized_counts.csv",
                           # clip_regular_file_out
                           "test_files/filter_annotate_binding_regions/clip_test1_normalized_counts_out.csv",
                           # results_list
                           [",x50_sample_1,x50_sample_2,x50_sample_3,x5_sample_1,x5_sample_2,x5_sample_3,gene_id,gene_type,gene_name,sub_gene_type,all_genes,all_sub_gene_types",
                            "chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693.2,TEC,4933401J01Rik,exon,,"],
                           # count_match
                           (2)),
                         # Test 4: Multiple matching events with multiple genes/sub-genes with a 
                         # normalized counts file, with regular CLIP files.
                         # in an exon region.
                         # sm_filtered_file
                         (("test_files/filter_annotate_binding_regions/sm_test4_normalized_counts.csv"),
                           # regular_gtf_file
                           ("test_files/filter_annotate_binding_regions/gtf_test4.csv"),
                           # sm_filtered_file_out
                           ("test_files/filter_annotate_binding_regions/sm_test2_normalized_counts_out.csv"),
                           # clip_regular_file
                           "test_files/filter_annotate_binding_regions/clip_test2_normalized_counts.csv",
                           # clip_regular_file_out
                           "test_files/filter_annotate_binding_regions/clip_test2_normalized_counts_out.csv",
                           # results_list
                           [",x50_sample_1,x50_sample_2,x50_sample_3,x5_sample_1,x5_sample_2,x5_sample_3,gene_id,gene_type,gene_name,sub_gene_type,all_genes,all_sub_gene_types",
                            "chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693.2,TEC,4933401J01Rik,exon,ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992,exon~3_UTR",
                            "chr1~100~300~-,100.0,100.3,100.2,50.1,50.2,50.4,,,,,,",
                            "chr1~3100~3300~+,100.0,100.3,100.2,50.1,50.2,50.4,,,,,,",
                            "chr2~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693,snRNA,gene,,,",
                            "chr6~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,,,,,,"],
                           # count_match
                           (6)),
                         ])

def test_filter_annotate_binding_regions(sm_filtered_file, regular_gtf_file, sm_filtered_file_out,
                                         clip_regular_file, clip_regular_file_out, results_list, count_match):
    """
        GIVEN a DeSeq2 produced file that was SM filtered, with or without additional
            regular, non-filtered CLIP DeSeq2 produced file, and a GTF annotation, the
            files are annotated and the regular CLIP file is filtered using the SM filtered
            file.
        WHEN the function annotates the DeSeq2 input binding regions with the GTF features,
            and filters the addition non-filtered CLIP DeSeq2 events (if input) and outputs
            them in a file.
        THEN the output file is checked for the correct output by
            counting the number of perfectly matching lines.
    """
    # Sets line_count to 0. This is the count of how many
    # lines match in the output and the expected output file.
    line_count = 0
    # Executes function.
    f_a_b_r.filter_annotate_binding_regions(sm_filtered_file, regular_gtf_file, sm_filtered_file_out,
                                            clip_regular_file, clip_regular_file_out)
    # Opens testing output files to iterate through to check against
    # the expected output in the the results_list.
    read1_opened = open(sm_filtered_file_out, 'r')
    # Checks each line of the output against the expected output
    # in the results_list. Counts should match.
    for line in read1_opened:
        # Cleans line to make it is easier to compare.
        line_clean = line.strip("\n")
        if line_clean in results_list:
            line_count += 1
    # Checks line counts.
    assert count_match == line_count
    read1_opened.close()