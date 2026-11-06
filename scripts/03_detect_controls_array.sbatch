#!/bin/bash
#SBATCH --job-name=str-detect-controls
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=96:00:00
# ^ The LAST alignment dominates and scales with the number of bases. A
#   control brought to roughly 18x needs about 70 hours at 16 cores on the
#   measured profile, so 24 hours kills every task mid alignment. Raise this
#   further, or lower the subsampling target, if your controls are deeper.
#SBATCH --partition=<your_partition>
#SBATCH --array=1-N%16
# ^ Replace N with the number of data lines in controls.tsv. The throttle is
#   what sets the wall clock: at 70 hours per sample, 110 controls take 80
#   days at %4 and 13 days at %24. Raise it as far as the disk allows, since
#   each concurrent task holds its own FASTQ and MAF.
# ^ Replace N with the number of lines in controls.tsv (excluding header):
#     tail -n +2 controls.tsv | wc -l
#   %4 limits to 4 concurrent jobs -- adjust based on resources available
#   on your partition (each task already uses 16 CPU / 64G).

# A job array: each task processes ONE line of the TSV (sample_id, bam_path,
# fastq_path). Required before `build-controls`, which reads
# results/controls/{sample_id}/{sample_id}.merged.vcf.

set -euo pipefail
mkdir -p logs

cd "${SLURM_SUBMIT_DIR}"

SAMPLES_TSV="controls.tsv"
# +1 to skip the TSV header line
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${SAMPLES_TSV}")
SAMPLE_ID=$(echo "${LINE}" | cut -f1)
BAM=$(echo "${LINE}" | cut -f2)
FASTQ=$(echo "${LINE}" | cut -f3)

echo "Control ${SLURM_ARRAY_TASK_ID}: ${SAMPLE_ID} (bam=${BAM}, fastq=${FASTQ})"

str-toolkit detect \
  --sample "${SAMPLE_ID}" \
  --bam "${BAM}" \
  --fastq "${FASTQ}" \
  --config config.yaml \
  --threads "${SLURM_CPUS_PER_TASK}" \
  -o results/controls/
