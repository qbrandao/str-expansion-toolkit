#!/bin/bash
#SBATCH --job-name=str-truthset
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --partition=<your_partition>

# Validates meiotic instability calls against the Platinum Pedigree tandem
# repeat truthset for CEPH 1463 (Kronenberg et al., Nat Methods 2025).
#
# The truthset is a set of TR genotypes retained because they segregate
# consistently with Mendelian inheritance. At those loci transmission is
# known to be consistent, so a parent-to-child difference called by the
# pipeline is most likely a false positive. This yields an empirical
# false-positive rate for the nearest-size transmitted-allele assignment.
#
# Download the truthset first (once):
#   aws s3 cp --no-sign-request \
#     s3://platinum-pedigree-data/variants/tr_truthset/GRCh38/ceph_1463_tandem_repeats.oa.vcf.gz \
#     references/
#
# Run AFTER 09_meiotic_instability.sbatch.
#
# IMPORTANT: run the inspect step first (see below) and confirm which INFO
# field carries the motif, rather than assuming it.

set -euo pipefail
mkdir -p logs

cd "${SLURM_SUBMIT_DIR}"

TRUTHSET="references/ceph_1463_tandem_repeats.oa.vcf.gz"

# Step 0, run once interactively rather than as a job:
#   str-toolkit validate-truthset --truthset "${TRUTHSET}" --inspect -o /dev/null

str-toolkit validate-truthset \
  --truthset "${TRUTHSET}" \
  --instability results/meiotic_instability.tsv \
  --window 25 \
  --tolerance 0 \
  -o results/meiotic_instability.truthset_annotated.tsv \
  --summary results/truthset_false_positive_rate.tsv
