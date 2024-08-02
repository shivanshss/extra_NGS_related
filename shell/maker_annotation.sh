#!/bin/bash

# README
# This script performs various bioinformatics analyses including BUSCO assessment, repeat modeling, 
# and repeat masking on a given genome. It also isolates and formats complex repeats for further 
# processing with MAKER. The script uses several tools such as BUSCO, RepeatModeler, RepeatMasker, 
# and custom Perl scripts for these tasks.

# Author: Bloodmark

# 1. Setup
# - Creates a working directory and navigates to it.
# - Defines input files for the genome, transcriptome, and proteins.
# - Runs BUSCO to assess the completeness of the genome assembly.

# 2. Repeat Modeling
# - Sets up environment variables for RepeatModeler.
# - Builds a database from the genome.
# - Runs RepeatModeler to identify repeat elements in the genome.

# 3. Repeat Masking
# - Sets up environment variables for RepeatMasker.
# - Runs RepeatMasker using a custom repeat library.
# - Runs RepeatMasker with species-specific repeat libraries.
# - Combines and processes repeat masker output files.

# 4. Complex Repeats Isolation and Formatting
# - Isolates complex repeats from the RepeatMasker output.
# - Reformats the isolated repeats for use with MAKER.

# Input details:
# - Genome file: /home/dagashe/shivanshss/new_genomes/tcas6_finalGapfillScaffCorr.fasta
# - Transcriptome file: /home/dagashe/shivanshss/download27/transcriptome_assembly/JNS_RNAseq/SOAP_assembly_fasta_files/SOAP_31.fasta
# - Proteins file: /home/dagashe/shivanshss/uniprot_sprot.fasta

# Output details:
# - BUSCO results: busco_tcas6 directory
# - RepeatModeler log: repeatmodeler.log
# - RepeatMasker output: repeatMasker-output directory
# - Combined RepeatMasker output: full_mask directory
# - Isolated and reformatted complex repeats: full_mask.out.complex.reformat.gff3

# Script Efficiency Suggestions:
# - Consider parallelizing the BUSCO and RepeatMasker steps using GNU Parallel or a similar tool.
# - Use variables for file paths and parameters to avoid redundancy and make the script easier to modify.
# - Clean up intermediate files as you go to save disk space.
# - Add error handling to catch and report issues at each step.

# Example of running the script:
# bash your_script_name.sh

# Note: Ensure all the necessary tools and dependencies are installed and accessible.

# Code starts below:


mkdir /home/dagashe/shivanshss/final_maker

cd /home/dagashe/shivanshss/final_maker

MY_GENOME=/home/dagashe/shivanshss/new_genomes/tcas6_finalGapfillScaffCorr.fasta
MY_TRANSCRIPTOME=/home/dagashe/shivanshss/download27/transcriptome_assembly/JNS_RNAseq/SOAP_assembly_fasta_files/SOAP_31.fasta
MY_PROTEINS=/home/dagashe/shivanshss/uniprot_sprot.fasta

/home/dagashe/shivanshss/tools/busco/bin/busco -i /home/dagashe/shivanshss/new_genomes/tcas6_finalGapfillScaffCorr.fasta -l insecta_odb10 -o busco_tcas6 --long -m genome --augustus --offline --config /home/dagashe/shivanshss/tools/busco/config/config_tcas6


export PATH=/home/dagashe/shivanshss/tools/RepeatModeler-2.0.3:$PATH
export PERL5LIB=/softwares/perl5.30.0/lib/
export PERL5LIB=/home/dagashe/shivanshss/perl5/lib
export PERL5LIB=/home/dagashe/shivanshss/perl5/lib/perl5

BuildDatabase -name genome_db -engine ncbi $MY_GENOME

RepeatModeler -pa 32 -engine ncbi -database genome_db 2>&1 | tee repeatmodeler.log

#################################################################################################################################################
REPEATMASKER_LIB_DIR=/home/dagashe/shivanshss/tools/RepeatMasker/Libraries

REPEATMASKER_MATRICES_DIR=/home/dagashe/shivanshss/tools/RepeatMasker/Matrices

export PATH=/home/dagashe/shivanshss/tools/RepeatMasker/:$PATH


mkdir repeatMasker-tcas_model

RepeatMasker -pa 24 -e ncbi -gccalc 
[-lib combined-conseni.fa.classified #concatenated consensi.fa.classified and RepBase files]
-dir repeatMasker-model $MY_GENOME

mkdir repeatMasker-output

repeat_species=arthropoda

RepeatMasker -pa 24 -e ncbi -gccalc -a -species $repeat_species -dir repeatMasker-output $MY_GENOME 

mkdir full_mask

gunzip repeatMasker-output/*.cat.gz 

cat repeatMasker-output/*.cat.gz > full_mask/full_mask.cat

cd full_mask

ProcessRepeats -species arthropoda full_mask.cat

# create GFF3
rmOutToGFF3.pl full_mask/full_mask.out > full_mask/full_mask.out.gff3


# isolate complex repeats
grep -v -e "Satellite" -e ")n" -e "-rich" full_mask.out.gff3 > full_mask.out.complex.gff3

# reformat to work with MAKER
cat full_mask.out.complex.gff3 | perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' > full_mask.out.complex.reformat.gff3

complex_repeats=full_mask.out.complex.reformat.gff3


