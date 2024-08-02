#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Script to Extract Contigs from a Scaffold FASTA File

Description:
This script reads a scaffold FASTA file and extracts contigs separated by 'N' characters. Each extracted contig is written to a new FASTA file with the format `>contig_<number>`.

Author: Bloodmark

Input Details:
- The input is a scaffold FASTA file specified by the `-i` option.
- The script expects the input file path as an argument.

Output Details:
- The output is a FASTA file containing the extracted contigs, specified by the `-o` option.
- Each contig is named `>contig_<number>`.

How to Use:
1. Save the script as `contig_from_scaffold.py`.
2. Run the script from the command line with the input and output file paths as arguments:
python contig_from_scaffold.py -i <input_scaffold_fasta> -o <output_contig_fasta>

3. The script will generate an output FASTA file with the extracted contigs.

Example:
python contig_from_scaffold.py -i scaffold.fasta -o contigs.fasta


This will read `scaffold.fasta` and create `contigs.fasta` with the extracted contigs.

Note:
- Ensure that the input file follows the expected FASTA format.
- The script separates contigs based on sequences of 'N' characters in the scaffold.

Script Functioning:
1. Parses command-line arguments to get the input and output file paths.
2. Reads the scaffold FASTA file and concatenates sequences.
3. Splits the concatenated sequence into contigs at 'N' characters.
4. Writes each contig to the output FASTA file with the format `>contig_<number>`.
"""

from Bio import SeqIO
import getopt
import sys
import re

def usage():
    print "Usage: python contig_from_scaffold.py -i <input_scaffold_fasta> -o <output_contig_fasta>"

try:
    options, remainder = getopt.getopt(sys.argv[1:], 'i:o:h')
except getopt.GetoptError as err:
    print str(err)
    usage()
    sys.exit(2)

input_file = None
output_file = None

for opt, arg in options:
    if opt in ('-i'):
        input_file = arg
    elif opt in ('-o'):
        output_file = arg
    elif opt in ('-h'):
        usage()
        sys.exit()

if not input_file or not output_file:
    usage()
    sys.exit(2)

out = open(output_file, 'w')

sequence = ''.join([str(record.seq).strip() for record in SeqIO.parse(input_file, "fasta")])

m = re.sub('[nN]+', '\n', sequence).split('\n')

for i in range(len(m)):
    if m[i]:  # Only write non-empty sequences
        out.write('>contig_' + str(i + 1) + '\n')
        out.write(m[i] + '\n')

out.close()

