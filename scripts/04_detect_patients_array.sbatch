#!/bin/bash
#SBATCH --job-name=str-detect-patients
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=<your_partition>
#SBATCH --array=1-N%4
# ^ Replace N with the number of lines in patients.tsv (excluding header):
#     tail -n +2 patients.tsv | wc -l

set -euo pipefail
mkdir -p logs

cd "${SLURM_SUBMIT_DIR}"

SAMPLES_TSV="patients.tsv"
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${SAMPLES_TSV}")
SAMPLE_ID=$(echo "${LINE}" | cut -f1)
BAM=$(echo "${LINE}" | cut -f2)
FASTQ=$(echo "${LINE}" | cut -f3)

echo "Patient ${SLURM_ARRAY_TASK_ID}: ${SAMPLE_ID} (bam=${BAM}, fastq=${FASTQ})"

str-toolkit detect \
  --sample "${SAMPLE_ID}" \
  --bam "${BAM}" \
  --fastq "${FASTQ}" \
  --config config.yaml \
  --threads "${SLURM_CPUS_PER_TASK}" \
  -o results/patients/
