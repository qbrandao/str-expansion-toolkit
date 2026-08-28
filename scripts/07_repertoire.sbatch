#!/bin/bash
#SBATCH --job-name=str-repertoire
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=<your_partition>

# Builds the genome-wide VNTR repertoire (paper Methods 2.8-2.9, Results 3.2):
# every locus seen in the control cohort, classified by genomic location
# (subtelomeric / paracentromeric / 5' region / exonic / intronic /
# intergenic_other) and by motif length category.
#
# Lightweight: only parses the merged VCFs already produced by
# 03_detect_controls_array.sbatch. Memory is set higher than compare because
# the whole repertoire is held as a dataframe before writing.
#
# Run AFTER all control detection array jobs have completed.

set -euo pipefail
mkdir -p logs

cd "${SLURM_SUBMIT_DIR}"
str-toolkit repertoire \
  --controls-dir results/controls/ \
  --samples-list controls.tsv \
  --genes-bed genes.bed.gz \
  --exons-bed MANE_Select_exons.bed.gz \
  --promoter-bp 2000 \
  -o results/repertoire.tsv \
  --summary results/repertoire_summary.tsv
