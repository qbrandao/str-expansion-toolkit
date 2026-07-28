#!/bin/bash
#SBATCH --job-name=str-compare
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --partition=<ta_partition>

# Léger : parse les VCF fusionnés des patients + results/controls.json,
# annote gène/feature, trie. À lancer APRÈS 04_detect_patients_array.sbatch
# ET 05_build_controls.sbatch (cf. submit_full_pipeline.sh).

set -euo pipefail
mkdir -p logs

cd "${SLURM_SUBMIT_DIR}"
str-toolkit compare \
  --patients-dir results/patients/ \
  --controls-json results/controls.json \
  --genes-bed genes.bed.gz \
  --exons-bed MANE_Select_exons.bed.gz \
  -x 0 \
  -o results/report.tsv
