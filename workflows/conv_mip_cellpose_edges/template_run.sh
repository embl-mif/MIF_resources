#!/bin/bash

#SBATCH --job-name=conv_mip_cellpose_edges
#SBATCH --partition=htc-el8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j_%u.out
#SBATCH --error=logs/%x_%j_%u.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=user.name@embl.es # <---------------- MODIFY

set -euo pipefail

# Load Nextflow.
module load Nextflow/25.10.1

# Nextflow launch/work directory. Create it before submitting this script.
cd "/scratch/$USER/conv_mip_cellpose_edges" # <---------- MODIFY IF NEEDED

# Exporting these variables makes them visible to Nextflow subprocesses.
export NF_PATH="/home/$USER/projects/mif-resources/workflows/conv_mip_cellpose_edges/main.nf" # <--- MODIFY
export WDIR="/PATH/TO/WORK/DIRECTORY" <----------------MODIFY TO YOUR DESIRED WORK LOCATION
export CSV_PATH="/PATH/TO/batch.csv" # <----------------- MODIFY
# export OUTDIR="/scratch/$USER/conv_mip_cellpose_edges/results" # <--- MODIFY

nextflow run "$NF_PATH" \
    -work-dir "$WDIR" \
    --input_csv "$CSV_PATH" \
    -resume
