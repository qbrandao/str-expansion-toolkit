#!/bin/bash
#SBATCH --job-name=str-detect-p01
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=<your_partition>

# Full detect (default tools: VAMOS + tandem-genotypes + LongTR, plus
# merge) for ONE patient. clair3 (VAMOS) and last-train/lastal on the
# whole genome are the slowest/most resource-intensive steps -- see the
# README for resource guidance.

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
