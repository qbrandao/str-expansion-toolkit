#!/bin/bash
#SBATCH --job-name=str-detect-ceph
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=<your_partition>
#SBATCH --array=1-4%2

# Runs `detect` on the four CEPH 1463 quartet members used for the germline
# instability analysis (paper Methods 2.7):
#   NA12877 (father), NA12878 (mother), NA12879 (daughter), NA12882 (son)
#
# Required before 09_meiotic_instability.sbatch, which reads
# results/ceph/{sample_id}/{sample_id}.merged.vcf.
#
# ceph.tsv must be a TSV with header: sample_id, bam_path, fastq_path
# %2 limits to 2 concurrent jobs -- raise if your partition allows.

set -euo pipefail
mkdir -p logs

cd "${SLURM_SUBMIT_DIR}"

SAMPLES_TSV="ceph.tsv"
# +1 to skip the TSV header line
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${SAMPLES_TSV}")
SAMPLE_ID=$(echo "${LINE}" | cut -f1)
BAM=$(echo "${LINE}" | cut -f2)
FASTQ=$(echo "${LINE}" | cut -f3)

echo "CEPH sample ${SLURM_ARRAY_TASK_ID}: ${SAMPLE_ID} (bam=${BAM}, fastq=${FASTQ})"

str-toolkit detect \
  --sample "${SAMPLE_ID}" \
  --bam "${BAM}" \
  --fastq "${FASTQ}" \
  --config config.yaml \
  --threads "${SLURM_CPUS_PER_TASK}" \
  -o results/ceph/
