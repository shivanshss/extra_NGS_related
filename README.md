# extra_NGS_related
This repository contains simple codes and functions often used in NGS (next Gen seq) and TGS (Third gen seq) analyses

## 1. Bulk download sequences from NCBI using python

folder called python_download_sequences

NCBI Sequence Downloader

This Python script is designed to download nucleotide sequences from the NCBI database using accession numbers. It uses the Biopython library to interact with the NCBI Entrez system and fetch sequences in FASTA format. Accession numbers and corresponding names are provided in a text file, and the script downloads the sequences and saves them as FASTA files.

Usage:

    Ensure you have Python 3.x installed on your system.
    Install the Biopython library if you haven't already. You can do this using the command: pip install biopython.
    Replace "your_email_address" in the script with your actual email address. This is required by NCBI Entrez.
    Prepare a text file with accession numbers and names in tab-separated format. Example: accession_list.txt.
    Set the input_file variable in the script to the path of your input file (e.g., input_file = "/path/to/accession_list.txt").
    Specify the output_directory where you want to save the downloaded sequences.
    Run the script using the command: python script_name.py.

Script Functionality:

    The script utilizes the provided Biopython library and NCBI Entrez API to fetch sequences by accession numbers.
    Accession numbers and corresponding names are read from the input file.
    The script attempts to download sequences for each accession number and saves them as individual FASTA files in the specified output directory.
    If any errors occur during the download process, appropriate messages are displayed.

Note:

    This script assumes that you have access to the NCBI Entrez system and that the provided accession numbers are valid.
    Make sure to replace "your_email_address" with your actual email address to comply with NCBI's usage policy.
    The script requires Biopython to be installed, which can be done using the command: pip install biopython.

Feel free to use this script to conveniently download nucleotide sequences from NCBI using accession numbers for your research purposes.

## 2. Split multi-sample multi-chromosome vcf file by chromosomes

vcf_to_chr.sh file

Split VCF by Chromosome

This script is designed to split an uncompressed VCF (Variant Call Format) file into separate files based on chromosomes. It employs tabix and parallel processing to efficiently extract entries for each chromosome present in the VCF.

Usage:

    Ensure that the script has execution permissions. If not, you can grant it using the command: chmod +x split_vcf_by_chr.sh.
    Execute the script by providing the uncompressed VCF file as an argument: ./split_vcf_by_chr.sh input.vcf.

Script Functionality:

    The script first checks if the provided file is an uncompressed VCF. If not, an error message is displayed.
    It extracts the list of chromosomes present in the VCF and stores them in a text file.
    The script then sorts the VCF entries based on chromosome and position.
    The sorted entries are compressed and indexed using bgzip and tabix commands, resulting in a compressed and indexed VCF file.
    The script then proceeds to extract separate VCF files for each chromosome using parallel processing.

Note:

    This script assumes that the input VCF file is uncompressed.
    The extracted VCF files for each chromosome are named as ${input.vcf}_${chromosome}.vcf.

Feel free to use this script to split your VCF files into chromosome-specific files for further analysis or downstream processing.
