#!/bin/bash

# Constants
README_CONTENT="### Few pointers for Project Organization\n- Use GitHub for collaborating, sharing, and version control of your scripts!\n- Add your raw data, backup data, and any large files >50 Mb to .gitignore file.\n- Use a pipeline management tool such as snakemake or nextflow\n- Document each step of the metagenome analysis process.\n- Store intermediate files for reproducibility."

# Function to create directory with validation
create_directory() {
    if [ ! -d "$1" ]; then
        mkdir "$1"
        echo "Created directory: $1"
    else
        echo "Directory '$1' already exists. Skipping..."
    fi
}

# User input for project name
read -p "Enter the name of your project (avoid spaces, use underscores): " project_name

# Create project directory
create_directory "$project_name"

# Create subdirectories
create_directory "$project_name/data"
create_directory "$project_name/data/raw_reads"
create_directory "$project_name/data/processed_reads"

create_directory "$project_name/assembly"
create_directory "$project_name/assembly/assembly_results"

create_directory "$project_name/functional_annotation"
create_directory "$project_name/functional_annotation/annotation_results"

create_directory "$project_name/metagenome_analysis"
create_directory "$project_name/metagenome_analysis/metabarcoding"
create_directory "$project_name/metagenome_analysis/metagenomic_analysis"

create_directory "$project_name/results"
create_directory "$project_name/results/analysis_results"

create_directory "$project_name/docs"
create_directory "$project_name/docs/protocols"
create_directory "$project_name/docs/analysis_notes"

create_directory "$project_name/scripts"
create_directory "$project_name/scripts/assembly_scripts"
create_directory "$project_name/scripts/annotation_scripts"
create_directory "$project_name/scripts/metagenome_analysis_scripts"

create_directory "$project_name/tmp"
create_directory "$project_name/tmp/assembly_tmp"
create_directory "$project_name/tmp/annotation_tmp"
create_directory "$project_name/tmp/metagenome_analysis_tmp"

# Create readme files
echo "$README_CONTENT" > "$project_name/README.md"
echo "This directory is for raw metagenome reads." > "$project_name/data/raw_reads/readme.md"
echo "This directory is for processed metagenome reads." > "$project_name/data/processed_reads/readme.md"

echo "This directory contains metagenome assembly results." > "$project_name/assembly/assembly_results/readme.md"

echo "This directory contains functional annotation results." > "$project_name/functional_annotation/annotation_results/readme.md"

echo "This directory contains metabarcoding analysis results." > "$project_name/metagenome_analysis/metabarcoding/readme.md"
echo "This directory contains metagenomic analysis results." > "$project_name/metagenome_analysis/metagenomic_analysis/readme.md"

echo "This directory contains analysis results for each step." > "$project_name/results/analysis_results/readme.md"

echo "This directory contains protocols for various steps." > "$project_name/docs/protocols/readme.md"
echo "This directory contains analysis notes for each step." > "$project_name/docs/analysis_notes/readme.md"

echo "This directory contains scripts for metagenome assembly." > "$project_name/scripts/assembly_scripts/readme.md"
echo "This directory contains scripts for functional annotation." > "$project_name/scripts/annotation_scripts/readme.md"
echo "This directory contains scripts for metagenome analysis." > "$project_name/scripts/metagenome_analysis_scripts/readme.md"

# Success message
echo "Project Initialization completed successfully!"
tree "$project_name"

