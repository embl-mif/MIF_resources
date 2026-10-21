#!/bin/bash

#SBATCH --job-name=endo_morph
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
cd "/scratch/$USER/endo_morph" # <---------- MODIFY IF NEEDED

# Exporting these variables makes them visible to Nextflow subprocesses.
export NF_PATH="/home/$USER/projects/mif_resources/workflows/endo_morph/main.nf" # <--- MODIFY
export LIB_PATH="/home/$USER/projects/mif_resources/lib" # <--- MODIFY
export WDIR="/PATH/TO/WORK/DIRECTORY" # <---------------- MODIFY TO YOUR DESIRED WORK LOCATION
export CSV_PATH="/PATH/TO/batch.csv" # <----------------- MODIFY
# Gaps (empty area) are set per row in the CSV, e.g. endomorph_empty_channels="Alexa_488 Alexa_647"
# (cell-body markers) and endomorph_empty_nucleus_channel=HOECHST_33342_clean (see README.md).
export MDL="/g/mif-hd/shared/cellpose_models/cpsam" # <--- MODIFY IF USING A DIFFERENT MODEL

# The merged cells/summary CSVs are saved next to the input csv.
export OUT_DIR="$(dirname "$CSV_PATH")/merged" # <--- MODIFY IF DESIRED

nextflow run "$NF_PATH" \
    -lib "$LIB_PATH" \
    -work-dir "$WDIR" \
    --input_csv "$CSV_PATH" \
    --merged_output_dir "$OUT_DIR" \
    --cellpose_model_path "$MDL" \
    -resume
