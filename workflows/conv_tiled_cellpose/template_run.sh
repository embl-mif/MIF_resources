#!/bin/bash

#SBATCH --job-name=complete_conv_tiled_cellpose
#SBATCH --partition=htc-el8
#SBATCH --mem=16G
#SBATCH --time=72:00:00
#SBATCH --output=logs/%x_%j_%u.out
#SBATCH --error=logs/%x_%j_%u.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=user.name@embl.es # <---------------- MODIFY

set -euo pipefail

# Load Nextflow.
module load Nextflow/25.10.1

# Nextflow launch directory. Create it before submitting this script.
cd "/scratch/$USER/complete_conv_tiled_cellpose" # <---------- MODIFY IF NEEDED

# Exporting these variables makes them visible to Nextflow subprocesses.
export NF_PATH="/home/$USER/projects/mif_resources/workflows/complete_conv_tiled_cellpose/main.nf" # <--- MODIFY
export LIB_PATH="/home/$USER/projects/mif_resources/lib" # <--- MODIFY
export WDIR="/PATH/TO/WORK/DIRECTORY" # <---------------- MODIFY TO YOUR DESIRED WORK LOCATION
export CSV_PATH="/PATH/TO/batch.csv" # <----------------- MODIFY
export MDL="/g/mif-hd/shared/cellpose_models/cpsam" # <--- MODIFY IF USING A DIFFERENT MODEL

# Optional pipeline-wide parameters. Uncomment the ones you want to override.
#   --outdir              default location for the converted zarr/mask (default: results)
#   --resolution_level    pyramid level for tiling/segmentation when a CSV row has none (default: 0)
#   --label_offset_step   label-ID offset step between tiles (default: 100000)
#   --multicut_beta       LABEL_SOLVER edge-cut trade-off (default: 0.5)
EXTRA_ARGS=(
    # --outdir "/scratch/$USER/complete_conv_tiled_cellpose/results"
    # --resolution_level 0
    # --label_offset_step 100000
    # --multicut_beta 0.5
)

nextflow run "$NF_PATH" \
    -lib "$LIB_PATH" \
    -work-dir "$WDIR" \
    --input_csv "$CSV_PATH" \
    --cellpose_model_path "$MDL" \
    "${EXTRA_ARGS[@]}" \
    -resume
