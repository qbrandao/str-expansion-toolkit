#!/bin/bash
#SBATCH --job-name=str-build-controls
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --partition=<ta_partition>

# Léger : ne fait que parser les VCF fusionnés déjà produits par
# 03_detect_controls_array.sbatch (results/controls/{id}/{id}.merged.vcf).
# À lancer APRÈS que tous les jobs array de détection contrôles soient
# terminés (cf. submit_full_pipeline.sh pour l'enchaînement automatique).

set -euo pipefail
mkdir -p logs

cd "${SLURM_SUBMIT_DIR}"
str-toolkit build-controls \
  --controls-dir results/controls/ \
  --samples-list controls.tsv \
  -o results/controls.json
