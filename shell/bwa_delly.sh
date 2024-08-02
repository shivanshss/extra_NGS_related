#!/bin/bash

# README
# This script performs a series of bioinformatics analyses including alignment, SNP calling,
# structural variant detection, and visualization using tools such as NUCmer, BWA, Samtools,
# Delly, BCFtools, Breakdancer, and GRIDSS. The script is divided into sections for different
# tasks, and the input files are referenced accordingly. The outputs are also described for each step.

# Author: Bloodmark

# 1. Show coordinates and SNPs from NUCmer delta file
# Input: /home/bloodmark/workarea/nucmer/NC_007417.3.delta
# Outputs:
#   - NC_007417.coords.txt: Coordinates of alignments
#   - NC_007417.snps.txt: SNPs detected

# 2. BWA Indexing and Alignment
# Inputs:
#   - /home/bloodmark/workarea/new_final_chr/tcas6_final_all.fasta
#   - /home/bloodmark/workarea/pubref/pubref_chromosomes/GCF_000002335.3_Tcas5.2_genomic.fna
# Outputs:
#   - alignment.sam: SAM file from BWA MEM alignment
#   - alignment.bam: BAM file converted from SAM
#   - alignment_sorted.bam: Sorted BAM file
#   - alignment_sorted.bam.bai: BAM index file

# 3. Structural Variant Detection with Delly
# Inputs:
#   - /home/bloodmark/workarea/pubref/pubref_chromosomes/LGY.fa
#   - lgy_alignment_sorted.bam
# Outputs:
#   - lgy_output.bcf: BCF file with structural variants
#   - lgy_output.vcf: VCF file with structural variants
#   - delly_output.bed: BED file with structural variants
#   - filtered_delly.bed: BED file with filtered structural variants (INS, DEL, DUP, INV, TRA)

# 4. Specific Chromosome Alignment and Structural Variant Detection
# Inputs:
#   - /home/bloodmark/workarea/new_final_chr/NC_007425.3_RagTag.fa (file1)
#   - /home/bloodmark/workarea/pubref/pubref_chromosomes/LG10.fa (file2)
# Outputs:
#   - ${name}_alignment.sam: SAM file from BWA MEM alignment
#   - ${name}_alignment.bam: BAM file converted from SAM
#   - ${name}_alignment_sorted.bam: Sorted BAM file
#   - ${name}_alignment_sorted.bam.bai: BAM index file
#   - ${name}_output.bcf: BCF file with structural variants
#   - ${name}_output.vcf: VCF file with structural variants
#   - ${name}_delly_output.bed: BED file with structural variants
#   - ${name}_filtered_delly.bed: BED file with filtered structural variants (INS, DEL, DUP, INV, TRA)

# 5. Run Breakdancer for Structural Variant Detection
# Inputs:
#   - ${name}_alignment_sorted.bam
# Outputs:
#   - ${name}_breakdancer.config: Breakdancer configuration file
#   - ${name}_breakdancer.sv: Breakdancer structural variants file
#   - ${name}_breakdancer.bed: BED file with structural variants

# 6. GRIDSS for Structural Variant Detection
# Inputs:
#   - ${name}_alignment_sorted.bam
# Outputs:
#   - ${name}_output.vcf.gz: Compressed VCF file with structural variants
#   - ${name}_structural_variants.bed: BED file with structural variants

# 7. Visualization with MUMmerplot
# Inputs: Various delta files
# Outputs: PNG images visualizing the alignments

# Script Efficiency Suggestions:
# - Consider parallelizing the BWA indexing and alignment steps for multiple files using GNU Parallel or a similar tool.
# - Use variables for file paths to avoid redundancy and make the script easier to modify.
# - Clean up intermediate files as you go to save disk space.
# - Combine consecutive commands where possible to reduce I/O operations.

# Example of running the script:
# bash your_script_name.sh

# Note: Uncomment the relevant sections to run the desired steps.

# Code starts below:



#show-coords -rcl /home/bloodmark/workarea/nucmer/NC_007417.3.delta > NC_007417.coords.txt

#show-snps -Clr /home/bloodmark/workarea/nucmer/NC_007417.3.delta > NC_007417.snps.txt

#bwa index /home/bloodmark/workarea/new_final_chr/tcas6_final_all.fasta
#bwa index /home/bloodmark/workarea/pubref/pubref_chromosomes/GCF_000002335.3_Tcas5.2_genomic.fna 

#bwa mem /home/bloodmark/workarea/new_final_chr/tcas6_final_all.fasta /home/bloodmark/workarea/pubref/pubref_chromosomes/GCF_000002335.3_Tcas5.2_genomic.fna > alignment.sam
#samtools view -Sb alignment.sam > alignment.bam
#samtools sort alignment.bam -o alignment_sorted.bam
#samtools index alignment_sorted.bam

#delly call -g /home/bloodmark/workarea/pubref/pubref_chromosomes/LGY.fa -o lgy_output.bcf lgy_alignment_sorted.bam

#bcftools view lgy_output.bcf > lgy_output.vcf

#bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\n' lgy_output.vcf > delly_output.bed
#grep -E 'INS|DEL|DUP|INV|TRA' delly_output.bed > filtered_delly.bed



file1=/home/bloodmark/workarea/new_final_chr/NC_007425.3_RagTag.fa
file2=/home/bloodmark/workarea/pubref/pubref_chromosomes/LG10.fa
name=chr10

#bwa index $file1
#bwa index $file2

#bwa mem $file1 $file2 > ${name}_alignment.sam
#samtools view -Sb ${name}_alignment.sam > ${name}_alignment.bam
#samtools sort ${name}_alignment.bam -o ${name}_alignment_sorted.bam
#samtools index ${name}_alignment_sorted.bam

#rm ${name}_alignment.bam ${name}_alignment.sam

#delly call -g $file1 -o ${name}_output.bcf ${name}_alignment_sorted.bam

#bcftools view ${name}_output.bcf > ${name}_output.vcf

#bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\n' ${name}_output.vcf > ${name}_delly_output.bed


#grep -E 'INS|DEL|DUP|INV|TRA' ${name}_delly_output.bed > ${name}_filtered_delly.bed

# Run Breakdancer
#breakdancer-max -d ${name}_breakdancer -g $file1 ${name}_alignment_sorted.bam

# Convert Breakdancer output to BED format
#perl /path/to/breakdancer/perl/bam2cfg.pl ${name}_alignment_sorted.bam > ${name}_breakdancer.config
#breakdancer-max -r 4 -t -q 10 ${name}_breakdancer.config > ${name}_breakdancer.sv
#perl /path/to/breakdancer/perl/sv2bed.pl ${name}_breakdancer.sv > ${name}_breakdancer.bed

#GRIDSS_JAR=/home/bloodmark/Tools/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar

#/home/bloodmark/Tools/gridss/gridss --reference $file1 --output ${name}_output.vcf.gz ${name}_alignment_sorted.bam --jar /home/bloodmark/Tools/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar --jvmheap "4g"

#bcftools query -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\n' ${name}_output.vcf.gz | awk 'BEGIN {OFS="\t"} {end = $2 + $4; print $1, $2, end, $3}' > ${name}_structural_variants.bed


mummerplot -p lgy -l lgy.delta --layout -t png --large --filter

mummerplot -p tcon -l tcas6_final_all_tcon.delta --layout -t png --large --filter

mummerplot -p lg2 -l NC_007417.3.delta --layout -t png --large --filter

mummerplot -p lg3 -l NC_007418.3.delta --layout -t png --large --filter

mummerplot -p lg4 -l NC_007419.2.delta --layout -t png --large --filter

mummerplot -p lg5 -l NC_007420.3.delta --layout -t png --large --filter

mummerplot -p lg6 -l NC_007421.3.delta --layout -t png --large --filter

mummerplot -p lg7 -l NC_007422.5.delta --layout -t png --large --filter

mummerplot -p lg8 -l NC_007423.3.delta --layout -t png --large --filter

mummerplot -p lg9 -l NC_007424.3.delta --layout -t png --large --filter

mummerplot -p lg10 -l NC_007425.3.delta --layout -t png --large --filter

mummerplot -p tfree -l tcas6_final_all_tfree.delta --layout -t png --large --filter



