"""test_remove_lower_regions.py- Unit tests for
    remove_lower_regions.py
"""
import pytest
# Imports all functions from remove_lower_regions.py.
import remove_lower_regions as r_b_r

"""
Tests get_sample_names() with different numbers of samples and missing data.

Data Types:
    sample_csv_file -- Text file containing the file name to sample names in
        the normalized counts file. The filename should not contain the directory
        structure. This is used to match the files to the samples in the normalized
        counts file. The file is in CSV format. No title required.
        Format- filename,deseq_sample_name,..   
    sample_names -- Dictionary of the sample names.
        {filename: [deseq_sample_name]}.
    message -- Message to be printed.       
    """

@pytest.mark.parametrize("sample_csv_file, sample_names, message", [
                         # Test 1: One line with three samples.
                         # sample_csv_file
                         (('test_files/remove_lower_regions/test1_sample_csv_file.csv'),
                           # sample_names
                           ({'name1.csv': ["s1","s2","s3"]}),
                           # message
                           ("")),
                         # Test 2: Multiple lines with three samples.
                         # sample_csv_file
                         (('test_files/remove_lower_regions/test2_sample_csv_file.csv'),
                           # sample_names
                           ({'name1.csv': ["s1","s2","s3"],
                             'name2.csv': ["s4","s5","s6"],
                             'name3.csv': ["s7","s8","s9"]}),
                           # message
                           ("")),         
                         # Test 3: Wrongly formatted file (tabs and samples
                         # not added).
                            # sample_csv_file
                            (('test_files/remove_lower_regions/test3_sample_csv_file.csv'),
                            # sample_names
                            ({'name1.csv\ts1\ts2\ts3': [],
                              'name2.csv': []}),
                            # message
                            ("ERROR: No sample names found in sample_csv_file or sample"
                             " names are not formatted correctly. System will continue"
                             " to run, but samples might be missing. Please check file.\n")),
                         # Test 4: Empty input file.
                            # sample_csv_file
                            (('test_files/remove_lower_regions/test40_sample_csv_file.csv'),
                            # sample_names
                            ({}),
                            # message
                            ("ERROR: No sample names found in sample_csv_file or sample"
                             " names are not formatted correctly. System will continue"
                             " to run, but samples might be missing. Please check file.\n")),
                         ])

def test_get_sample_names(sample_csv_file, sample_names, message, capsys):
    """
        GIVEN a sample CSV file for identify samples in the normalized 
            counts file.
        WHEN the function parses the file and adds the file name and the 
            corresponding sample names to a dictionary.
        THEN the output dictionary is checked for the key/value pairs.
    """
    assert r_b_r.get_sample_names(sample_csv_file) == (sample_names)
    # Checks the print statement.
    captured = capsys.readouterr()
    assert captured.out == message

"""
Tests get_average_counts() with different numbers input samples
    and binding regions.

Data Types:
    normalized_counts -- DESeq2 produced normalized counts file. 
        This file should contain all the normalized counts for all 
        protein concentrations.
    sample_names -- Dictionary of the sample names.
        {filename: [deseq_sample_name]}.   
    average_counts -- Dictionary of the average normalized counts 
        for each sample.
        {sample_name: {binding_region: average_normalized_count}}.
"""

@pytest.mark.parametrize("normalized_counts, sample_names, average_counts", [
                          # Test 1: One line with three samples.
                          # normalized_counts
                          (('test_files/remove_lower_regions/test1_normalized_counts.csv'),
                            # sample_names
                            ({'name1.csv': ["s1","s2","s3"]}),
                            # average_counts
                            ({'name1.csv': {'region_1': 2.0}})),
                          # Test 2: Multiple lines with nine samples.
                          # normalized_counts
                            (('test_files/remove_lower_regions/test2_normalized_counts.csv'),
                            # sample_names
                            ({'name1.csv': ["s1","s2","s3"],
                              'name2.csv': ["s5","s8","s10"],
                              'name3.csv': ["s6","s7","s9"]}),
                            # average_counts
                            ({'name1.csv': {'region_1': 2.0,
                                            'region_2': 6.0,
                                            'region_3': 1.0},
                              'name2.csv': {'region_1': 1.0,
                                            'region_2': 3.0,
                                            'region_3': 3.0},
                              'name3.csv': {'region_1': 2.0,
                                            'region_2': 4.0,
                                            'region_3': 6.0}})),
                        ])         

def test_get_average_counts(normalized_counts, sample_names, average_counts):
    """
        GIVEN a normalized counts file and a dictionary of sample names.
        WHEN the function parses the file and calculates the average
            normalized counts for each sample.
        THEN the output dictionary is checked for the key/value pairs.
    """
    assert r_b_r.get_average_counts(normalized_counts, sample_names) == (average_counts)

"""
Tests remove_lower_regions() with different numbers of binding regions.

Data Types:
    deseq_file_list -- DESeq2 produced file type. Should contain all files for 
        all protein concentrations. All the files shoule be the SM filtered
        and filtered with the control sample (i.e., 
        x50_vs_uni_significant_filt_annot.csv). The files should be listed 
        in the order desired to be filtered. This is is usually the lower 
        concentration to higher concentration. The file is in CSV format.

    lower_regions -- Dictionary of binding regions that were found in the
        lower concentration files.
        {"file_name": {binding_region: line}}. 
    all_regions -- Dictionary of all binding regions found in the files.
        {"file_name": {binding_region: line}}.
    title_dictionary -- Dictionary of the title lines for each file.
        Used to output the title lines to the output files.
        {"file_name": title_line}.
"""

@pytest.mark.parametrize("deseq_file_list, lower_regions, all_regions, title_dictionary", [
                          # Test 1: One identical line with three sample files.
                          # deseq_file_list
                          ((['test_files/remove_lower_regions/test1_deseq_file.csv',
                             'test_files/remove_lower_regions/test12_deseq_file.csv',
                             'test_files/remove_lower_regions/test13_deseq_file.csv']),
                            # lower_regions
                           ({'test1_deseq_file.csv': {'region_1': 'region_1,line1'},
                             'test12_deseq_file.csv': {},
                             'test13_deseq_file.csv': {}}),
                            # all_regions
                           ({'test1_deseq_file.csv': {'region_1': 'region_1,line1'},
                             'test12_deseq_file.csv': {'region_1': 'region_1,line1'},
                             'test13_deseq_file.csv': {'region_1': 'region_1,line1'}}),
                            # title_dictionary
                           ({'test1_deseq_file.csv': ',title_line',
                             'test12_deseq_file.csv': ',title_line',
                             'test13_deseq_file.csv': ',title_line'})),
                          # Test 2: Multiple lines with nine samples.
                          # deseq_file_list
                          ((['test_files/remove_lower_regions/test2_deseq_file.csv',
                             'test_files/remove_lower_regions/test22_deseq_file.csv',
                             'test_files/remove_lower_regions/test23_deseq_file.csv']),
                            # lower_regions
                            ({'test2_deseq_file.csv': {'region_1': 'region_1,line1',
                                                       'region_2': 'region_2,line2'},
                              'test22_deseq_file.csv': {'region_3': 'region_3,line3'},
                              'test23_deseq_file.csv': {'region_4': 'region_4,line4',
                                                        'region_5': 'region_5,line5'}}),
                            # all_regions
                            ({'test2_deseq_file.csv': {'region_1': 'region_1,line1',
                                                        'region_2': 'region_2,line2'},
                              'test22_deseq_file.csv': {'region_1': 'region_1,line1',
                                                        'region_2': 'region_2,line2',
                                                        'region_3': 'region_3,line3'},
                              'test23_deseq_file.csv': {'region_1': 'region_1,line1',
                                                        'region_2': 'region_2,line2',
                                                        'region_3': 'region_3,line3',
                                                        'region_4': 'region_4,line4',
                                                        'region_5': 'region_5,line5'}}),
                            # title_dictionary
                            ({'test2_deseq_file.csv': ',title_line',
                             'test22_deseq_file.csv': ',title_line',
                             'test23_deseq_file.csv': ',title_line'})),
                            ])

def test_remove_lower_regions(deseq_file_list, lower_regions, all_regions, title_dictionary):
    """
        GIVEN a list of DESeq2 files, a dictionary of binding regions
            from the lower concentration files, a dictionary of all
            binding regions, and a dictionary of title lines.
        WHEN the function parses the files and removes the binding
            regions from the lower concentration files.
        THEN the output dictionaries are checked for the key/value pairs.
    """
    assert r_b_r.remove_lower_regions(deseq_file_list) == (lower_regions, all_regions, 
                                                           title_dictionary)

"""
Tests output_file() with different numbers of binding regions and samples.

Data Types:
    average_counts -- Dictionary of the average normalized counts 
        for each sample.
        {sample_name: {binding_region: average_normalized_count}}
    lower_region -- Dictionary of binding regions that were found in the
        lower concentration files.
        {"file": {binding_region: line}} .
    all_regions -- Dictionary of all binding regions found in the files.
        {"file": {binding_region: line}}.
    title_dictionary -- Dictionary of the title lines for each file.
        Used to output the title lines to the output files.
        {"file": title_line}.
    output -- Output to add to the output file name.
    region_output_check_list -- List of CSV files of expected 
        nucleotide and motif secondary statistics.
    count_match_list -- List of number of lines that match 
        between the output and the first check file.    
"""

@pytest.mark.parametrize("average_counts, lower_regions, all_regions, \
                         title_dictionary, output, region_output_check_list, \
                         count_match_list", [
                         # Test 1: One sample line.
                         # average_counts
                         ({'name1.csv': {'region_1': 2.0}},
                          # lower_regions
                          ({'name1.csv': {'region_1': 'region_1,line1'}}),
                          # all_regions
                          ({'name1.csv': {'region_1': 'region_1,line1'}}),
                          # title_dictionary
                          ({'name1.csv': ',title_line'}),
                          # output
                          ('test_files/remove_lower_regions'),
                          # region_output_check_list
                          (['test_files/remove_lower_regions/name1_no_lower_check.csv',
                            'test_files/remove_lower_regions/name1_counts_check.csv']),
                          # count_match_list
                          ([2, 2])),
                          # Test 2: Multiple samples and lines. Also tests output formatting.
                         # average_counts
                         ({'name12.csv': {'region_1': 2.0,
                                          'region_2': 6.0,
                                          'region_3': 1.0},
                           'name22.csv': {'region_1': 1.0,
                                          'region_2': 3.0,
                                          'region_3': 3.0},
                            'name32.csv': {'region_1': 2.0,
                                           'region_2': 4.0,
                                           'region_3': 5.0}},
                          # lower_regions (below not realistic but used for testing)
                          ({'name12.csv': {'region_1': 'region_1,line1',
                                           'region_2': 'region_2,line2',
                                           'region_3': 'region_3,line3'},
                            'name22.csv': {'region_1': 'region_1,line1',
                                           'region_2': 'region_2,line2',
                                           'region_3': 'region_3,line3'},
                            'name32.csv': {'region_1': 'region_1,line1',
                                           'region_2': 'region_2,line2'}}),
                          # all_regions
                          ({'name12.csv': {'region_1': 'region_1,line1',
                                           'region_2': 'region_2,line2',
                                           'region_3': 'region_3,line3'},
                            'name22.csv': {'region_1': 'region_1,line1',
                                           'region_2': 'region_2,line2',
                                           'region_3': 'region_3,line3'},
                            'name32.csv': {'region_1': 'region_1,line1',
                                           'region_2': 'region_2,line2',
                                           'region_3': 'region_2,line3'}}),
                          # title_dictionary
                          ({'name12.csv': ',title_line',
                             'name22.csv': ',title_line',
                             'name32.csv': ',title_line'}),
                          # output
                          ('test_files/remove_lower_regions/'),
                          # region_output_check_list
                          (['test_files/remove_lower_regions/name12_no_lower_check.csv',
                            'test_files/remove_lower_regions/name12_counts_check.csv',
                            'test_files/remove_lower_regions/name22_no_lower_check.csv',
                            'test_files/remove_lower_regions/name22_counts_check.csv',
                            'test_files/remove_lower_regions/name32_no_lower_check.csv',
                            'test_files/remove_lower_regions/name32_counts_check.csv']),
                          # count_match_list
                          ([4, 4, 4, 4, 3, 3])), 
                         ])

def test_output_file(average_counts, lower_regions, all_regions, title_dictionary, output,
                      region_output_check_list, count_match_list):
    """
    GIVEN a dictionary of the average normalized counts for each sample,
        a dictionary of binding regions from the lower concentration files,
        a dictionary of all binding regions, a dictionary of title lines,
        and a output will be used to output the count information.
    WHEN the function organizes and outputs the information to the two
        output files.
    THEN the output files are checked for the correct information.
    """
    r_b_r.output_file(average_counts, lower_regions, all_regions, title_dictionary, output)
    # Counts used to check output.
    line_count = 0
    # Dictionary used to save files needed to be checked.
    file_dictionary = {}
    # Grabs file name and adds them to the file_list.
    for file in average_counts:
        # Saves the filename without the directory structure.
        file_name = file.split("/")[-1].replace(".csv", "")
        # Opens the right output files so content can be checked.
        if output[-1] == "/":
            deseq_file_name = (output + file_name + "_no_lower.csv")
            normalized_counts_name = (output + file_name + "_counts.csv")
        else:
            deseq_file_name = (output + "/" + file_name + "_no_lower.csv")
            normalized_counts_name = (output + "/" + file_name + "_counts.csv")
        # Saves files by file name.
        file_dictionary[file_name] = [deseq_file_name, normalized_counts_name]  
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

"""Tests process_regions() with different numbers of samples and binding regions.

Data Types:
    deseq_file_list -- DESeq2 produced file type. Should contain all files for 
        all protein concentrations. All the files shoule be the SM filtered
        and filtered with the control sample (i.e., 
        x50_vs_uni_significant_filt_annot.csv). The files should be listed 
        in the order desired to be filtered. This is is usually the lower 
        concentration to higher concentration. The file is in CSV format.
    normalized_counts -- DESeq2 produced normalized counts file. This file
        should contain all the normalized counts for all protein concentrations.
    sample_csv_file -- Text file containing the file name to sample names in
        the normalized counts file. The filename should not contain the directory
        structure. This is used to match the files to the samples in the normalized
        counts file. The file is in CSV format. No title required.
        Format- filename,deseq_sample_name,..
    output -- Output to add to the output file name.
    region_output_check_list -- List of CSV files of expected 
        nucleotide and motif secondary statistics.
    count_match_list -- List of number of lines that match 
        between the output and the first check file.
"""

@pytest.mark.parametrize("deseq_file_list, normalized_counts, sample_csv_file, output, \
                          region_output_check_list, count_match_list", [
                            # Test 1: One sample line.
                            # deseq_file_list
                         ((['test_files/remove_lower_regions/test3_deseq_file.csv',
                            'test_files/remove_lower_regions/test32_deseq_file.csv',
                            'test_files/remove_lower_regions/test33_deseq_file.csv']),
                         # normalized_counts
                         ('test_files/remove_lower_regions/test31_normalized_counts.csv'),
                         # sample_csv_file
                         ('test_files/remove_lower_regions/test31_sample_csv_file.csv'),
                         # output
                         ('test_files/remove_lower_regions/'),
                         # region_output_check_list
                         (['test_files/remove_lower_regions/test3_deseq_file_no_lower_check.csv',
                            'test_files/remove_lower_regions/test3_deseq_file_counts_check.csv',
                            'test_files/remove_lower_regions/test32_deseq_file_no_lower_check.csv',
                            'test_files/remove_lower_regions/test32_deseq_file_counts_check.csv',
                            'test_files/remove_lower_regions/test33_deseq_file_no_lower_check.csv',
                            'test_files/remove_lower_regions/test33_deseq_file_counts_check.csv']),
                         # count_match_list
                         ([2, 2, 1, 2, 1, 2])),
                         # Test 2: Multiple samples and lines. Also tests output formatting.
                         # deseq_file_list
                         ((['test_files/remove_lower_regions/test4_deseq_file.csv',
                            'test_files/remove_lower_regions/test42_deseq_file.csv',
                            'test_files/remove_lower_regions/test43_deseq_file.csv']),
                         # normalized_counts
                         ('test_files/remove_lower_regions/test4_normalized_counts.csv'),
                         # sample_csv_file
                         ('test_files/remove_lower_regions/test4_sample_csv_file.csv'),
                         # output
                         ('test_files/remove_lower_regions/'),
                         # region_output_check_list
                         (['test_files/remove_lower_regions/test4_deseq_file_no_lower_check.csv',
                            'test_files/remove_lower_regions/test4_deseq_file_counts_check.csv',
                            'test_files/remove_lower_regions/test42_deseq_file_no_lower_check.csv',
                            'test_files/remove_lower_regions/test42_deseq_file_counts_check.csv',
                            'test_files/remove_lower_regions/test43_deseq_file_no_lower_check.csv',
                            'test_files/remove_lower_regions/test43_deseq_file_counts_check.csv']),
                         # count_match_list
                         ([3, 3, 2, 4, 3, 6])),
                        ])

def test_process_regions(deseq_file_list, normalized_counts, sample_csv_file, output,
                         region_output_check_list, count_match_list):
     """
     GIVEN a list of DESeq2 files, a normalized counts file, a sample CSV file,
          and a output will be used to output the count information.
     WHEN the function organizes and outputs the information to the two
          output files.
     THEN the output files are checked for the correct information.
     """
     r_b_r.process_regions(deseq_file_list, normalized_counts, sample_csv_file, output)
     # Counts used to check output.
     line_count = 0
     # Dictionary used to save files needed to be checked.
     file_dictionary = {}
     # Grabs file name and adds them to the file_list.
     for file in deseq_file_list:
        # Saves the filename without the directory structure.
        file_name = file.split("/")[-1].replace(".csv", "")
        # Opens the right output files so content can be checked.
        if output[-1] == "/":
            deseq_file_name = (output + file_name + "_no_lower.csv")
            normalized_counts_name = (output + file_name + "_counts.csv")
        else:
            deseq_file_name = (output + "/" + file_name + "_no_lower.csv")
            normalized_counts_name = (output + "/" + file_name + "_counts.csv")
        # Saves files by file name.
        file_dictionary[file_name] = [deseq_file_name, normalized_counts_name]
     # Count used to open correct check file.
     count = 0
     # Loops through and checks output of each file. Four loops because the number of
     # iterations is short.
     for file_name in file_dictionary:
          file_list = file_dictionary[file_name]
          for file in file_list:
                # Counts used to check output.
                line_count = 0
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