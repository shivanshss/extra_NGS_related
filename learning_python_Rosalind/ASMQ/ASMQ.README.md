# DNA assembly Statistics Calculator

This script calculates N50, N75, L50, and total genome size from a collection of FASTA files. It reads a list of FASTA file names, processes the sequences, and outputs the statistics in a tabular format.

## Usage:

    Create an input file "fasta_files.txt" with a list of FASTA filenames, one per line.
    Change the script to and add input and output file names. (Changed to be interactive, it would ask for file names upon running)
    The output file will contain a table with FASTA file names, N50, N75, L50, and total genome size.

## Instructions:

    Save FASTA filenames in "fasta_files.txt".
    Run script: python script_name.py.
    Enter input filename (e.g., "fasta_files.txt") and output filename (e.g., "output_stats.txt").
    Check the output file for the statistics table.

## Example Input:


fasta_files.txt:
sequence1.fasta
sequence2.fasta

## Example Output:


output_stats.txt:
File Name     N50     N75     L50     Total Genome Size
sequence1.fasta 1500    800     2       6000
sequence2.fasta 1200    600     3       4800

## Note:

    Each input FASTA file should contain DNA sequences, following the standard FASTA format.
