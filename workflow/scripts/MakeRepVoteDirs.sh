#!/usr/bin/env bash

# 1. Create the main project directory
mkdir -p RepVote
cd RepVote

# 2. Create the core workflow hierarchy
mkdir -p workflow/rules
mkdir -p workflow/scripts

# 3. Create configuration and reference data folders
mkdir -p config
mkdir -p resources
mkdir -p .test
mkdir -p benchmarks

# 4. Initialize empty files for the workflow
touch workflow/Snakefile
touch workflow/rules/kmermaking.smk
touch workflow/rules/download.smk
touch workflow/rules/quantify.smk

# 5. Initialize configuration and documentation files
touch config/config.yaml
touch config/samples.tsv
touch README.md
touch LICENSE

# 6. (Optional) Create a .gitignore to keep the repo clean
cat <<EOF > .gitignore
.snakemake/
results/
data/
benchmarks/*
*.sif
.DS_Store
EOF

echo "RepVote structure created successfully."
