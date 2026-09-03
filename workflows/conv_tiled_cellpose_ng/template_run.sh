#!/bin/bash

#SBATCH --job-name=conv_tiled_cellpose_ng
#SBATCH --partition=htc-el8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j_%u.out
#SBATCH --error=logs/%x_%j_%u.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=user.name@embl.es # <---------------- MODIFY

set -euo pipefail

module load Nextflow/25.10.1

# Directory from which Nextflow will be launched.
cd "/scratch/$USER" # <----------- MODIFY IF NEEDED

export NF_PATH="/home/$USER/projects/mif_resources/workflows/conv_tiled_cellpose_ng/main.nf" # <--- MODIFY
export WDIR="/scratch/$USER/work" # <--- MODIFY
export CSV_PATH="/PATH/TO/batch.csv" # <---------------- MODIFY
export OUTDIR="/scratch/$USER/results" # <--- MODIFY

nextflow run "$NF_PATH" \
    -work-dir "$WDIR" \
    --input_csv "$CSV_PATH" \
    --outdir "$OUTDIR" \
    -resume
