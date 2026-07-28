#!/bin/bash
#SBATCH --job-name=str-detect-p01
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=<ta_partition>

# Detect complet (VAMOS + TRGT + tandem-genotypes + fusion) pour UN patient.
# clair3 (VAMOS) + last-train/lastal genome entier sont les étapes les plus
# lentes/gourmandes -> voir README pour le détail des ressources.

set -euo pipefail
mkdir -p logs

SAMPLE_ID="patient01"
BAM="/store/.../${SAMPLE_ID}.sorted.bam"
FASTQ="/store/.../${SAMPLE_ID}_Guppy_4.0.11_prom.merged.fastq.gz"

cd "${SLURM_SUBMIT_DIR}"
str-toolkit detect \
  --sample "${SAMPLE_ID}" \
  --bam "${BAM}" \
  --fastq "${FASTQ}" \
  --config config.yaml \
  --threads "${SLURM_CPUS_PER_TASK}" \
  -o results/patients/
