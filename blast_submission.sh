#!/bin/bash

# Bash script used to loop through a directory of FASTA files.
# Runs BLAST on each file and outputs to a file. All FASTA
# files should contain a ".fasta" extension.

# Input command line arguments.
###
# $1 -- Input directory to process all files in.
# $2 -- Taxonomic ID to BLAST against.
###
# Checks for "/" at end of file input. If not, adds
# "/" to end of file input.
if [[ $1 != */ ]]
then
    $1=$1"/"
fi
files=$1/*
# Loops through all files in input directory, runs BLAST,
# and outputs to file.
for file in $files
do
    # First checks if file is a ".fasta" or ".fa" file.
    if [[ $file != *.fasta ]] && [[ $file != *.fa ]]
    then
        continue
    fi  
    # Removes ".fasta" or ".fa" from file name.
    file=${file%.fa}
    file=${file%.fasta}
    output=$file"_blast.tsv"
    echo "Processing $file file..."
    # BLAST command.
    blastn -db nt -taxids $2 -query $file -outfmt 6 -out $output
    # Prints entire BLAST command
    echo "blastn -db nt -taxids $2 -query $file -outfmt 6 -out $output"
done
