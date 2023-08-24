#!/bin/bash

# Constants
README_CONTENT="### Few pointers for Project Organization\n- Use GitHub for collaborating, sharing, and version control of your scripts!\n- Add your raw data, backup data, and any large files >50 Mb to .gitignore file.\n- Use a pipeline management tool such as snakemake or nextflow\n- Document each step of the NGS analysis process.\n- Store intermediate files for reproducibility."

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
create_directory "$project_name/data/raw_data"
create_directory "$project_name/data/processed_data"

create_directory "$project_name/analysis"
create_directory "$project_name/analysis/preprocessing"
create_directory "$project_name/analysis/alignment"
create_directory "$project_name/analysis/variant_calling"
create_directory "$project_name/analysis/post_processing"

create_directory "$project_name/results"
create_directory "$project_name/results/analysis_results"

create_directory "$project_name/docs"
create_directory "$project_name/docs/protocols"
create_directory "$project_name/docs/analysis_notes"

create_directory "$project_name/scripts"
create_directory "$project_name/scripts/preprocessing_scripts"
create_directory "$project_name/scripts/alignment_scripts"
create_directory "$project_name/scripts/variant_calling_scripts"
create_directory "$project_name/scripts/post_processing_scripts"

create_directory "$project_name/tmp"
create_directory "$project_name/tmp/preprocessing_tmp"
create_directory "$project_name/tmp/alignment_tmp"
create_directory "$project_name/tmp/variant_calling_tmp"
create_directory "$project_name/tmp/post_processing_tmp"

# Create readme files
echo "$README_CONTENT" > "$project_name/README.md"
echo "This directory is for raw NGS data." > "$project_name/data/raw_data/readme.md"
echo "This directory is for processed NGS data." > "$project_name/data/processed_data/readme.md"

echo "This directory contains preprocessing analysis results." > "$project_name/analysis/preprocessing/readme.md"
echo "This directory contains alignment analysis results." > "$project_name/analysis/alignment/readme.md"
echo "This directory contains variant calling analysis results." > "$project_name/analysis/variant_calling/readme.md"
echo "This directory contains post-processing analysis results." > "$project_name/analysis/post_processing/readme.md"

echo "This directory contains analysis results for each step." > "$project_name/results/analysis_results/readme.md"

echo "This directory contains protocols for various steps." > "$project_name/docs/protocols/readme.md"
echo "This directory contains analysis notes for each step." > "$project_name/docs/analysis_notes/readme.md"

echo "This directory contains scripts for preprocessing." > "$project_name/scripts/preprocessing_scripts/readme.md"
echo "This directory contains scripts for alignment." > "$project_name/scripts/alignment_scripts/readme.md"
echo "This directory contains scripts for variant calling." > "$project_name/scripts/variant_calling_scripts/readme.md"
echo "This directory contains scripts for post-processing." > "$project_name/scripts/post_processing_scripts/readme.md"

# Success message
echo "Project Initialization completed successfully!"
tree "$project_name"

