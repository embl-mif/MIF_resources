#!/bin/bash

#SBATCH --job-name=nf_launcher
#SBATCH --partition=htc-el8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j_%u.out
#SBATCH --error=logs/%x_%j_%u.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=user.name@embl.es # <----------------MODIFY

# Load nextflow
module load Nextflow/25.10.1

# Navigate to repo
cd /home/$USER/projects # <----------------MODIFY TO YOUR LOCAL REPO PATH

#Using 'export' ensures these are visible to Nextflow sub-processes
export CSV_PATH="/PATH/TO/batch.csv"           # <----------------MODIFY

# Run Nextflow
nextflow run mif_resources/workflows/PYMIF2ZARR/main.nf \
    --input_csv  "$CSV_PATH" \
    -resume

