#!/bin/bash
#SBATCH --job-name=str-newbam-longtr
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --partition=<your_partition>

# STEP 1 of 2 for a newly received BAM: the fast check.
#
# Runs LongTR alone on the BAM as it stands. Nothing is realigned, clair3
# does not run, and LAST does not run, so this finishes in a fraction of the
# time of the full pipeline and still exercises the parts most likely to be
# misconfigured: the reference, the region catalog, the contig naming, and
# the depth.
#
# Run scripts/inspect_bam.sh first and read its sections 1, 2, 4 and 6.
#
# If this job produces a VCF with a reasonable number of genotyped loci,
# move on to scripts/13_test_new_bam_full.sbatch.

set -euo pipefail
mkdir -p logs

SAMPLE_ID="newbam01"
BAM="/CHANGE/ME/sample.bam"

cd "${SLURM_SUBMIT_DIR}"

str-toolkit detect \
  --sample "${SAMPLE_ID}" \
  --bam "${BAM}" \
  --config config.yaml \
  --tools longtr \
  --threads "${SLURM_CPUS_PER_TASK}" \
  -o results/new_bam_test/

# --- what to look at once the job finishes ----------------------------
OUT="results/new_bam_test/${SAMPLE_ID}"
VCF="${OUT}/${SAMPLE_ID}.longtr.vcf.gz"

echo
echo "=== LongTR output summary ==="
if [[ ! -f "${VCF}" ]]; then
  echo "no VCF produced, read the .err file above" >&2
  exit 1
fi

N_TOTAL="$(zcat "${VCF}" | grep -vc '^#' || true)"
N_GT="$(zcat "${VCF}" | awk '!/^#/ && $10 !~ /^\.\/\.|^\.\|\./ {n++} END {print n+0}')"
echo "records in the VCF        : ${N_TOTAL}"
echo "records with a genotype   : ${N_GT}"
if [[ "${N_TOTAL}" -gt 0 ]]; then
  echo "genotyped fraction        : $(awk -v a="${N_GT}" -v b="${N_TOTAL}" 'BEGIN{printf "%.1f %%", 100*a/b}')"
fi
echo
echo "A genotyped fraction near zero with a non-empty catalog almost always"
echo "means the contig naming of the BAM and of the regions BED disagree, or"
echo "that the depth is below LongTR --min-reads. Check section 2 and section 5"
echo "of the inspect_bam.sh report."
echo
echo "=== merged VCF ==="
ls -l "${OUT}"/*.merged.vcf 2>/dev/null || echo "(no merged VCF, expected with a single tool)"
echo
echo "=== first 3 genotyped records ==="
zcat "${VCF}" | awk '!/^#/ && $10 !~ /^\.\/\.|^\.\|\./' | head -3 | cut -c1-240
