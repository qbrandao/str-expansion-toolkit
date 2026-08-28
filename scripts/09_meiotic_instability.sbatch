#!/bin/bash
#SBATCH --job-name=str-meiotic
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=<your_partition>

# Germline (meiotic) instability from the four CEPH 1463 parent-offspring
# duos (paper Methods 2.10, Results 3.3), in a 2x2 design crossing parent sex
# with offspring sex.
#
# Sex chromosomes are excluded by default: nearest-size allele matching
# assumes two comparable alleles per locus in both duo members, which fails
# where either is hemizygous (NA12882 on X and Y).
#
# Two outputs:
#   - meiotic_instability.tsv        one row per (duo, locus, tool)
#   - meiotic_instability_by_duo.tsv one row per (duo, tool, location, motif)
# The PER-DUO file is the correct unit of replication for comparing duo
# types: loci within a duo are not independent observations (shared genetic
# background and sequencing run), so comparing duo types from the per-locus
# file would be pseudo-replication.
#
# Run AFTER 08_detect_ceph_array.sbatch has fully completed.

set -euo pipefail
mkdir -p logs

cd "${SLURM_SUBMIT_DIR}"
str-toolkit meiotic-instability \
  --duos duos_ceph.tsv \
  --data-dir results/ceph/ \
  --genes-bed genes.bed.gz \
  --exons-bed MANE_Select_exons.bed.gz \
  --promoter-bp 2000 \
  -o results/meiotic_instability.tsv \
  --summary results/meiotic_instability_by_duo.tsv
