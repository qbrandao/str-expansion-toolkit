#!/bin/bash
#SBATCH --job-name=str-build-controls
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --partition=<your_partition>

# Lightweight: only parses the merged VCFs already produced by
# 03_detect_controls_array.sbatch (results/controls/{id}/{id}.merged.vcf).
# Run AFTER all control detection array jobs have completed (see
# submit_full_pipeline.sh for automatic chaining).

set -euo pipefail
mkdir -p logs

cd "${SLURM_SUBMIT_DIR}"
str-toolkit build-controls \
  --controls-dir results/controls/ \
  --samples-list controls.tsv \
  -o results/controls.json
