"""test_filter_annotate_binding_regions_bed.py- Unit tests for
    filter_annotate_binding_regions_bed.py
"""
import pytest
# Imports all functions from join_binding_regions.py.
import filter_annotate_binding_regions_bed as f_a_b_r_b

"""
Tests add_events_to_dictionary() with different BED and GTF files.

Data Types:
    bed_file -- BED file with binding regions. This file
        would be produced from a separate experiment where there
        where no SM Input samples.
    regular_gtf_file -- Gencode/Ensembl/UCSC provided GTF file. Should be the GTF
        annotation for the target alignment sequence. The file is used to annotate
        the binding regions.
    region_dictionary -- Dictionary containing the event name (chromosome information)
        as the key and the entire BED line with gene information as the value.
        {"chromosome_gtf,cordinate1_gtf,cordinate2_gtf,strand_gtf": 
        "line_clean=gene_line=sub_gene_line"}      
    """

@pytest.mark.parametrize("bed_file, regular_gtf_file, region_dictionary", [
                         # Test 1: One matching event with BED file
                         # in an exon region.
                         # bed_file
                         (('test_files/filter_annotate_binding_regions_bed/test1.bed'),
                           # regular_gtf_file
                           ('test_files/filter_annotate_binding_regions_bed/gtf_test1.csv'),
                          # region_dictionary
                          ({'chr1~100~300~+': 'chr1,100,300,X,15,+,0.001=ENSMUSG00000102693.2~TEC~4933401J01Rik=exon'})),
                         # Test 2: One matching event with BED file
                         # in an exon region. Contains multiple additional lines in
                         # BED file and counts file that do not match.
                         # bed_file
                         (('test_files/filter_annotate_binding_regions_bed/test2.bed'),
                           # regular_gtf_file
                           ('test_files/filter_annotate_binding_regions_bed/gtf_test2.csv'),
                          # region_dictionary
                          ({'chr1~100~300~+': 'chr1,100,300,X,15,+,0.001=ENSMUSG00000102693.2~TEC~4933401J01Rik=exon',
                            'chr1~100~300~-': 'chr1,100,300,X,15,-,0.001==',
                            'chr1~3100~3300~+': 'chr1,3100,3300,X,15,+,0.001==',
                            'chr2~100~300~+': 'chr2,100,300,X,15,+,0.001=='})),
                         # Test 3: Multiple matching event with BED file
                         # in multiple, single annotated regions. Contains other additional lines in
                         # BED file and counts file that do not match.
                         # bed_file
                         (('test_files/filter_annotate_binding_regions_bed/test2.bed'),
                           # regular_gtf_file
                           ('test_files/filter_annotate_binding_regions_bed/gtf_test3.csv'),
                          # region_dictionary
                          ({'chr1~100~300~+': 'chr1,100,300,X,15,+,0.001=ENSMUSG00000102693.2~TEC~4933401J01Rik=exon',
                            'chr1~100~300~-': 'chr1,100,300,X,15,-,0.001=ENSMUSG00000089699.2~lncRNA~Gm1992=',
                            'chr1~3100~3300~+': 'chr1,3100,3300,X,15,+,0.001==',
                            'chr2~100~300~+': 'chr2,100,300,X,15,+,0.001=ENSMUSG00000102693~snRNA~gene=utr'})),
                         # Test 4: Multiple matching event with BED file
                         # in multiple, multiple overlapping annotated regions. Some annotated regions
                         # are smaller and/or larger than binding region. Contains multiple 
                         # additional lines in BED file and counts file that do not match.
                         # bed_file
                         (('test_files/filter_annotate_binding_regions_bed/test2.bed'),
                           # regular_gtf_file
                           ('test_files/filter_annotate_binding_regions_bed/gtf_test4.csv'),
                          # region_dictionary
                          ({'chr1~100~300~+': 'chr1,100,300,X,15,+,0.001=ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992=exon~3utr~utr',
                            'chr1~100~300~-': 'chr1,100,300,X,15,-,0.001==',
                            'chr1~3100~3300~+': 'chr1,3100,3300,X,15,+,0.001==',
                            'chr2~100~300~+': 'chr2,100,300,X,15,+,0.001=ENSMUSG00000102693~snRNA~gene=utr'}))
                         ])

def test_add_events_to_dictionary(bed_file, regular_gtf_file, region_dictionary):
    """
        GIVEN a BED file with binding regions that needs to be annotated.
        WHEN the function finds overlapping genes and sub-gene types that
            overlap with the binding region and matches. The information
            is then added to a dictionary.
        THEN the output dictionary is checked for the correct region and 
            gene information. The title line string from the input file 
            is also checked for correctness.
    """
    assert f_a_b_r_b.add_events_to_dictionary(bed_file, regular_gtf_file) == (region_dictionary)

"""
Tests annotate_output_file() with different combinations of BED files and GTF annotations.

Data Types:
    bed_file -- BED file with binding regions. This file
        would be produced from a separate experiment where there
        where no SM Input samples.
    regular_gtf_file -- Gencode/Ensembl/UCSC provided GTF file. Should be the GTF
        annotation for the target alignment sequence. The file is used to annotate
        the binding regions.
    bed_file_out -- GTF Annotated BED file with binding regions. This file
        would be produced from a separate experiment where there
        where no SM Input samples.
    results_list -- List of results expected in function produced files.
    count_match -- Expected count of matching lines in output and premade
        testing files                 
    """

@pytest.mark.parametrize("bed_file, regular_gtf_file, bed_file_out, \
                          results_list, count_match", [
                         # Test 1: One matching event with a normalized counts file, no regular CLIP files.
                         # in an exon region.
                         # bed_file
                         (("test_files/filter_annotate_binding_regions_bed/test1.bed"),
                           # regular_gtf_file
                           ("test_files/filter_annotate_binding_regions_bed/gtf_test1.csv"),
                           # bed_file_out
                           ("test_files/filter_annotate_binding_regions_bed/test1_out.bed"),
                           # results_list
                           ["chr1,100,300,X,15,+,0.001,ENSMUSG00000102693.2,TEC,4933401J01Rik,exon,,"],
                           # count_match
                           (1)),
                         # Test 2: Multiple matching events with multiple genes/sub-genes with a 
                         # normalized counts file, no regular CLIP files.
                         # in an exon region.
                         # bed_file
                         (("test_files/filter_annotate_binding_regions_bed/test2.bed"),
                           # regular_gtf_file
                           ("test_files/filter_annotate_binding_regions_bed/gtf_test3.csv"),
                           # bed_file_out
                           ("test_files/filter_annotate_binding_regions_bed/test2_out.bed"),
                           # results_list
                           ['chr1,100,300,X,15,+,0.001,ENSMUSG00000102693.2,TEC,4933401J01Rik,exon,,',
                            'chr1,100,300,X,15,-,0.001,ENSMUSG00000089699.2,lncRNA,Gm1992,,,',
                            'chr1,3100,3300,X,15,+,0.001,,,,,,',
                            'chr2,100,300,X,15,+,0.001,ENSMUSG00000102693,snRNA,gene,utr,,'],
                           # count_match
                           (4)),
                         # Test 3: Multiple matching events with multiple genes/sub-genes with a 
                         # normalized counts file, no regular CLIP files.
                         # in an exon region.
                         # bed_file
                         (("test_files/filter_annotate_binding_regions_bed/test2.bed"),
                           # regular_gtf_file
                           ("test_files/filter_annotate_binding_regions_bed/gtf_test4.csv"),
                           # bed_file_out
                           ("test_files/filter_annotate_binding_regions_bed/test3_out.bed"),
                           # results_list
                           ["chr1,100,300,X,15,+,0.001,ENSMUSG00000102693.2,TEC,4933401J01Rik,exon,ENSMUSG00000102693.2~TEC~4933401J01Rik~ENSMUSG00000089699.2~lncRNA~Gm1992,exon~3utr~utr",
                            "chr1,100,300,X,15,-,0.001,,,,,,",
                            "chr1,3100,3300,X,15,+,0.001,,,,,,",
                            "chr2,100,300,X,15,+,0.001,ENSMUSG00000102693,snRNA,gene,utr,,"],
                           # count_match
                           (4)),
                         ])

def test_annotate_output_file(bed_file, regular_gtf_file, bed_file_out, \
                              results_list, count_match):
    """
        GIVEN a BED file that needs to be annotated using a GTF file and output.
        WHEN the function annotates the BED file with the GTF features and outputs
            them in a file.
        THEN the output file is checked for the correct output by
            counting the number of perfectly matching lines.
    """
    # Sets line_count to 0. This is the count of how many
    # lines match in the output and the expected output file.
    line_count = 0
    # Executes function.
    f_a_b_r_b.annotate_output_file(bed_file, regular_gtf_file, bed_file_out)
    # Opens testing output files to iterate through to check against
    # the expected output in the the results_list.
    read1_opened = open(bed_file_out, 'r')
    # Checks each line of the output against the expected output
    # in the results_list. Counts should match.
    for line in read1_opened:
        # Cleans line to make it is easier to compare.
        line_clean = line.strip("\n").replace("\t",",")
        if line_clean in results_list:
            line_count += 1
    # Checks line counts.
    assert count_match == line_count