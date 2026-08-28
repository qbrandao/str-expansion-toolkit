#!/bin/bash
#SBATCH --job-name=str-somatic
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=<your_partition>

# Somatic (mitotic) instability from per-read heterogeneity within single
# samples (paper Methods 2.11, Results 3.4).
#
# Needs no family structure: every sample already run through `detect` can
# contribute -- controls, patients and CEPH individuals alike. Build
# all_samples.tsv accordingly (a TSV with at least a sample_id column;
# see the note below on directory layout).
#
# Reads the tools' RAW per-sample outputs ({id}.longtr.vcf.gz and
# {id}.tandem_genotypes.tsv), NOT the merged VCF, since per-read detail is
# collapsed during merging. VAMOS is not used here: it reports per-haplotype
# consensus rather than individual-read measurements.
#
# NOTE ON --detect-dir: this expects ONE directory containing all
# {sample_id}/ subdirectories. If controls, patients and CEPH samples were
# written to separate output directories (results/controls/, results/patients/,
# results/ceph/), either run this script once per cohort with a matching
# samples list, or symlink all sample directories into a single directory
# first. Memory is set higher than the other analysis steps because per-read
# distributions across a whole cohort are held in memory before writing.
#
# Run AFTER the relevant `detect` jobs have completed.

set -euo pipefail
mkdir -p logs

cd "${SLURM_SUBMIT_DIR}"
str-toolkit somatic-instability \
  --samples-list all_samples.tsv \
  --detect-dir results/all/ \
  --genes-bed genes.bed.gz \
  --exons-bed MANE_Select_exons.bed.gz \
  --promoter-bp 2000 \
  --min-off-allele-reads 3 \
  -o results/somatic_instability.tsv
