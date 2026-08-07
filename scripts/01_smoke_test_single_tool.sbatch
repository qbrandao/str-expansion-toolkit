#!/bin/bash
#SBATCH --job-name=str-smoke-longtr
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=<your_partition>

# Quick end-to-end test on a SINGLE tool (LongTR: minimap2 + LongTR
# genotyping) to validate that config.yaml, reference paths, and
# micromamba environments are correct, before running the full pipeline
# (default tool set, much slower).

set -euo pipefail
mkdir -p logs

SAMPLE_ID="patient01"
FASTQ="/store/.../${SAMPLE_ID}_Guppy_4.0.11_prom.merged.fastq.gz"

cd "${SLURM_SUBMIT_DIR}"
str-toolkit detect \
  --sample "${SAMPLE_ID}" \
  --fastq "${FASTQ}" \
  --config config.yaml \
  --tools longtr \
  --threads "${SLURM_CPUS_PER_TASK}" \
  -o results/smoke_test/
