#!/bin/bash
#SBATCH --job-name=str-smoke-trgt
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=<ta_partition>

# Test rapide de bout en bout sur UN SEUL outil (TRGT, le plus rapide des 3 :
# minimap2 + trgt genotype) pour valider que config.yaml, les chemins de
# référence et les environnements micromamba sont corrects, avant de lancer
# le pipeline complet (3 outils, beaucoup plus long).

set -euo pipefail
mkdir -p logs

SAMPLE_ID="patient01"
FASTQ="/store/.../${SAMPLE_ID}_Guppy_4.0.11_prom.merged.fastq.gz"

cd "${SLURM_SUBMIT_DIR}"
str-toolkit detect \
  --sample "${SAMPLE_ID}" \
  --fastq "${FASTQ}" \
  --config config.yaml \
  --tools trgt \
  --threads "${SLURM_CPUS_PER_TASK}" \
  -o results/smoke_test/
