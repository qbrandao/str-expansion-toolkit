#!/bin/bash
# Downloads the Adotto catalog (hg38) and converts it to the format
# expected by LongTR: chrom, start (1-based), end, motif[,motif2], name.
#
# Usage: ./scripts/prepare_longtr_bed.sh /mnt/references
set -euo pipefail

OUTDIR="${1:-/mnt/references}"
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# LongTR/README.md #tr-region-bed-file requires gawk (3-argument match(),
# a GNU extension not supported by mawk/classic POSIX awk).
if ! awk --version 2>/dev/null | grep -qi "gnu awk"; then
  echo "ERROR: gawk is required (the default awk does not appear to be GNU Awk)." >&2
  echo "       Try 'module load gawk' or 'micromamba install gawk'."              >&2
  exit 1
fi

RAW_BED="adotto_repeats.hg38.bed"
CLEAN_BED="longtr_adotto_hg38_clean.bed"

if [ ! -f "${RAW_BED}" ]; then
  wget -q https://zenodo.org/records/7987365/files/adotto_repeats.hg38.bed.gz
  gunzip adotto_repeats.hg38.bed.gz
fi

# Adotto catalog columns: chrom, start (0-based BED), end, INFO
# (INFO = "ID=xxx;MOTIFS=motif1,motif2;STRUC=..."). Converted in a single
# pass: start -> 1-based, motif and ID extraction.
gawk -F'\t' 'BEGIN {OFS="\t"}
  {
    start = $2 + 1                      # BED 0-based -> LongTR 1-based
    info  = $4

    if (match(info, /MOTIFS=([^;]+)/, m)) {
      motif = m[1]
    } else {
      motif = info                      # fallback: column already a raw motif
    }

    if (match(info, /ID=([^;]+)/, idm)) {
      name = idm[1]
    } else {
      name = info
    }

    print $1, start, $3, motif, name
  }
' "${RAW_BED}" > "${CLEAN_BED}"

echo "OK: ${OUTDIR}/${CLEAN_BED} ($(wc -l < "${CLEAN_BED}") loci)"
echo "Set this path in config.yaml under longtr.regions_bed"
