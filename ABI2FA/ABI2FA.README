# ABI Sequencing Consensus Generator

This script generates combined consensus FASTA sequences from pairs of ABI sequencing files.

## Introduction

This Python script takes pairs of ABI sequencing files as input and generates combined consensus sequences in FASTA format for each pair. It uses the Biopython library for sequence parsing and alignment.
Prerequisites

    Python 3.x
    Biopython library (pip install biopython)

## Usage

    Create a text file (input_list.txt) containing pairs of ABI file paths in each row, separated by a space.

    Example input_list.txt:

path/to/abi1/file1.ab1 path/to/abi1/file2.ab1
path/to/abi2/file1.ab1 path/to/abi2/file2.ab1

## Run the script:

python generate_consensus.py

    The script will generate combined consensus FASTA sequences for each pair of ABI files and save them in separate directories.

## Input

The script expects a text file containing pairs of ABI file paths, with one pair per row. Each row should contain two file paths separated by a space.

## Output

For each row in the input file, the script creates a directory named row_X (where X is the row number) within the specified output directory. Inside each row directory, a file named consensus.fasta will contain the combined consensus sequence in FASTA format.

## Example

Suppose you have the following input files:

path/to/abi1/file1.ab1 path/to/abi1/file2.ab1
path/to/abi2/file1.ab1 path/to/abi2/file2.ab1

Running the script will generate the following structure in the output directory:

output_consensus/
|-- row_1/
|   `-- consensus.fasta
|-- row_2/
|   `-- consensus.fasta

The consensus.fasta files will contain the combined consensus sequences for each row.
