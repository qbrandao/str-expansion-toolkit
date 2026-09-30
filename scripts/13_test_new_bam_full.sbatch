#!/bin/bash
#SBATCH --job-name=str-newbam-full
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=36:00:00
#SBATCH --partition=<your_partition>

# STEP 2 of 2 for a newly received BAM: the full three-tool run and the merge.
#
# Only worth launching once scripts/12_test_new_bam_longtr.sbatch has produced
# a sensible VCF, since the two slow stages here are hard to debug:
#   clair3 on the whole genome (VAMOS branch), several hours
#   last-train then lastal on the whole genome (tandem-genotypes branch),
#   the longest stage by far, which is why the time limit is 36 hours
#
# tandem-genotypes needs reads rather than an alignment. With no --fastq, the
# reads are extracted from the BAM once (primary alignments only, methylation
# tags kept) and reused. Supply --fastq instead whenever the original FASTQ is
# available: extraction from an aligned BAM loses any read that failed to
# align, and LAST would otherwise never see them.
#
# Every stage is idempotent: resubmitting after a timeout picks up where the
# job stopped, so raise --time rather than starting over. The flip side is
# that an output produced under a wrong configuration is silently reused, so
# delete the sample directory after changing config.yaml.

set -euo pipefail
mkdir -p logs

SAMPLE_ID="newbam01"
BAM="/CHANGE/ME/sample.bam"
# Leave FASTQ empty to extract the reads from the BAM instead.
FASTQ=""

cd "${SLURM_SUBMIT_DIR}"

ARGS=(
  --sample "${SAMPLE_ID}"
  --bam "${BAM}"
  --config config.yaml
  --threads "${SLURM_CPUS_PER_TASK}"
  -o results/new_bam_test/
)
[[ -n "${FASTQ}" ]] && ARGS+=(--fastq "${FASTQ}")

# Add --realign here if inspect_bam.sh section 2 showed another reference
# build, or another aligner or preset. That needs a FASTQ.

str-toolkit detect "${ARGS[@]}"

# --- what to look at once the job finishes ----------------------------
OUT="results/new_bam_test/${SAMPLE_ID}"
MERGED="${OUT}/${SAMPLE_ID}.merged.vcf"

echo
echo "=== per-tool outputs ==="
for f in "${OUT}/${SAMPLE_ID}_assembly.hap1.vcf" \
         "${OUT}/${SAMPLE_ID}_assembly.hap2.vcf" \
         "${OUT}/${SAMPLE_ID}.tandem_genotypes.tsv" \
         "${OUT}/${SAMPLE_ID}.longtr.vcf.gz"; do
  if [[ -f "${f}" ]]; then
    if [[ "${f}" == *.gz ]]; then n="$(zcat "${f}" | grep -vc '^#' || true)"
    elif [[ "${f}" == *.vcf ]]; then n="$(grep -vc '^#' "${f}" || true)"
    else n="$(grep -vc '^#' "${f}" || true)"; fi
    printf '  %-52s %8s records\n' "$(basename "${f}")" "${n}"
  else
    printf '  %-52s %8s\n' "$(basename "${f}")" "MISSING"
  fi
done

echo
echo "=== merged VCF ==="
if [[ ! -f "${MERGED}" ]]; then
  echo "no merged VCF produced, read the .err file above" >&2
  exit 1
fi
echo "records: $(grep -vc '^#' "${MERGED}")"
echo
echo "loci by number of supporting tools:"
# INFO/SOURCES lists one entry per haplotype or allele (vamos_hap1,
# longtr_allele2, ...), so the tool family is the entry with its trailing
# _hapN / _alleleN removed, and the same tool must not be counted twice.
grep -v '^#' "${MERGED}" | grep -o 'SOURCES=[^;]*' | sed 's/SOURCES=//' | awk -F',' '
{
  delete seen; n = 0
  for (i = 1; i <= NF; i++) {
    fam = $i; sub(/_[^_]*$/, "", fam)
    if (!(fam in seen)) { seen[fam] = 1; n++ }
  }
  count[n]++
}
END { for (k = 1; k <= 4; k++) if (k in count) printf "  %d tool(s): %7d loci\n", k, count[k] }'
echo
echo "A large majority of single-tool loci is expected: the catalogs of the"
echo "three tools do not cover the same repeats. What matters for the merge is"
echo "that the 2-tool and 3-tool counts are not zero, which would mean the"
echo "coordinate reconciliation matched nothing at all."
echo
echo "=== first 3 merged records ==="
grep -v '^#' "${MERGED}" | head -3 | cut -c1-240
