"""test_filter_annotate_binding_regions.py- Unit tests for
    filter_annotate_binding_regions.py
"""
import pytest
# Imports all functions from join_binding_regions.py.
import filter_annotate_binding_regions as f_a_b_r

"""
Tests add_events_to_dictionary() with different DeSeq2 annotation files and GTF files.

Data Types:
    filtered_file -- DeSeq2 produced SM counts file (CSV) for a single protein
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
                         # DeSeq2 file and counts file that do not match.
                         # filtered_file
                         (('test_files/filter_annotate_binding_regions/sm_test2_normalized_counts.csv'),
                           # regular_gtf_file
                           ('test_files/filter_annotate_binding_regions/gtf_test3.csv'),
                          # region_dictionary
                          ({'chr1~100~300~+': 'chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik=exon',
                            'chr1~100~300~-': 'chr1~100~300~-,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000089699.2~lncRNA~Gm1992=',
                            'chr1~3100~3300~+': 'chr1~3100~3300~+,100.0,100.3,100.2,50.1,50.2,50.4==',
                            'chr2~100~300~+': 'chr2~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693~snRNA~gene=utr'}),
                          # title_line
                          ',x50_sample_1,x50_sample_2,x50_sample_3,x5_sample_1,x5_sample_2,x5_sample_3'),
                         # Test 4: Multiple matching event with a normalized counts file
                         # in multiple, multiple overlapping annotated regions. Some annotated regions
                         # are smaller and/or larger than binding region. Contains multiple 
                         # additional lines in DeSeq2 file and counts file that do not match.
                         # filtered_file
                         (('test_files/filter_annotate_binding_regions/sm_test2_normalized_counts.csv'),
                           # regular_gtf_file
                           ('test_files/filter_annotate_binding_regions/gtf_test4.csv'),
                          # region_dictionary
                          ({'chr1~100~300~+': 'chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=exon~3utr~utr',
                            'chr1~100~300~-': 'chr1~100~300~-,100.0,100.3,100.2,50.1,50.2,50.4==',
                            'chr1~3100~3300~+': 'chr1~3100~3300~+,100.0,100.3,100.2,50.1,50.2,50.4==',
                            'chr2~100~300~+': 'chr2~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4=ENSMUSG00000102693~snRNA~gene=utr'}),
                          # title_line
                          ',x50_sample_1,x50_sample_2,x50_sample_3,x5_sample_1,x5_sample_2,x5_sample_3'),
                         # Test 5: Same test as four but with DeSeq2 differential expression file, instead
                         # of normalized counts. Multiple matching event with a normalized counts file
                         # in multiple, multiple overlapping annotated regions. Some annotated regions
                         # are smaller and/or larger than binding region. Contains multiple 
                         # additional lines in DeSeq2 file and counts file that do not match.
                         # filtered_file
                         (('test_files/filter_annotate_binding_regions/sm_test3_differential_expression.csv'),
                           # regular_gtf_file
                           ('test_files/filter_annotate_binding_regions/gtf_test4.csv'),
                          # region_dictionary
                          ({'chr1~100~300~+': 'chr1~100~300~+,200,5,1,4,0.005,0.05=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=exon~3utr~utr',
                            'chr1~100~300~-': 'chr1~100~300~-,200,5,1,4,0.005,0.05==',
                            'chr1~3100~3300~+': 'chr1~3100~3300~+,200,5,1,4,0.005,0.05==',
                            'chr2~100~300~+': 'chr2~100~300~+,200,5,1,4,0.005,0.05=ENSMUSG00000102693~snRNA~gene=utr'}),
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
    results_list -- List of results expected in function produced files.
    count_match -- Expected count of matching lines in output and premade
        testing files                 
    """

@pytest.mark.parametrize("sm_filtered_file, regular_gtf_file, sm_filtered_file_out, sample_keyword, \
                          clip_regular_file, clip_regular_file_out, results_list, count_match", [
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
                           # sample_keyword
                           "SM",
                           # clip_regular_file
                           None,
                           # clip_regular_file_out
                           None,
                           # results_list
                           [",x50_sample_1,x50_sample_2,x50_sample_3,x5_sample_1,x5_sample_2,x5_sample_3,gene_id,gene_type,gene_name,sub_gene_type,all_genes,all_sub_gene_types",
                            "chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693.2,TEC,4933401J01Rik,exon,ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992,exon~3utr~utr",
                            "chr1~100~300~-,100.0,100.3,100.2,50.1,50.2,50.4,,,,,,",
                            "chr1~3100~3300~+,100.0,100.3,100.2,50.1,50.2,50.4,,,,,,",
                            "chr2~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693,snRNA,gene,utr,,"],
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
                           # sample_keyword
                           "CLIP",
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
                           # sample_keyword
                           "CLIP",
                           # clip_regular_file
                           "test_files/filter_annotate_binding_regions/clip_test2_normalized_counts.csv",
                           # clip_regular_file_out
                           "test_files/filter_annotate_binding_regions/clip_test2_normalized_counts_out.csv",
                           # results_list
                           [",x50_sample_1,x50_sample_2,x50_sample_3,x5_sample_1,x5_sample_2,x5_sample_3,gene_id,gene_type,gene_name,sub_gene_type,all_genes,all_sub_gene_types",
                            "chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693.2,TEC,4933401J01Rik,exon,ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992,exon~3utr~utr",
                            "chr1~100~300~-,100.0,100.3,100.2,50.1,50.2,50.4,,,,,,",
                            "chr1~3100~3300~+,100.0,100.3,100.2,50.1,50.2,50.4,,,,,,",
                            "chr2~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693,snRNA,gene,utr,,"],
                           # count_match
                           (5)),
                         ])

def test_output_file(sm_filtered_file, regular_gtf_file, sm_filtered_file_out, sample_keyword,
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
    f_a_b_r.output_file(sm_filtered_file, regular_gtf_file, sm_filtered_file_out, sample_keyword,
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
                            "chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693.2,TEC,4933401J01Rik,exon,ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992,exon~3utr~utr",
                            "chr1~100~300~-,100.0,100.3,100.2,50.1,50.2,50.4,,,,,,",
                            "chr1~3100~3300~+,100.0,100.3,100.2,50.1,50.2,50.4,,,,,,",
                            "chr2~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693,snRNA,gene,utr,,"],
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
                            "chr1~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693.2,TEC,4933401J01Rik,exon,ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992,exon~3utr~utr",
                            "chr1~100~300~-,100.0,100.3,100.2,50.1,50.2,50.4,,,,,,",
                            "chr1~3100~3300~+,100.0,100.3,100.2,50.1,50.2,50.4,,,,,,",
                            "chr2~100~300~+,100.0,100.3,100.2,50.1,50.2,50.4,ENSMUSG00000102693,snRNA,gene,utr,,"],
                           # count_match
                           (5)),
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