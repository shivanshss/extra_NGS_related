# extra_NGS_related
This repository contains simple codes and functions often used in NGS (next Gen seq) and TGS (Third gen seq) analyses

## 1. Bulk download sequences from NCBI using python

folder called python_download_sequences

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
