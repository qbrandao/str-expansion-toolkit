#!/bin/bash
# Télécharge le catalogue Adotto (hg38) et le convertit au format attendu
# par LongTR : chrom, start(1-based), end, motif[,motif2], nom.
#
# Usage : ./scripts/prepare_longtr_bed.sh /mnt/references
set -euo pipefail

OUTDIR="${1:-/mnt/references}"
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# LongTR/README.md #tr-region-bed-file exige gawk (match() à 3 arguments,
# extension GNU non supportée par mawk/awk POSIX classique).
if ! awk --version 2>/dev/null | grep -qi "gnu awk"; then
  echo "ERREUR : gawk requis (le awk par défaut ne semble pas être GNU Awk)." >&2
  echo "         Essaie 'module load gawk' ou 'micromamba install gawk'."     >&2
  exit 1
fi

RAW_BED="adotto_repeats.hg38.bed"
CLEAN_BED="longtr_adotto_hg38_clean.bed"

if [ ! -f "${RAW_BED}" ]; then
  wget -q https://zenodo.org/records/7987365/files/adotto_repeats.hg38.bed.gz
  gunzip adotto_repeats.hg38.bed.gz
fi

# Colonnes du catalogue Adotto : chrom, start(0-based BED), end, INFO
# (INFO = "ID=xxx;MOTIFS=motif1,motif2;STRUC=..."). On convertit directement
# en un seul passage : start -> 1-based, extraction du motif et de l'ID.
gawk -F'\t' 'BEGIN {OFS="\t"}
  {
    start = $2 + 1                      # BED 0-based -> LongTR 1-based
    info  = $4

    if (match(info, /MOTIFS=([^;]+)/, m)) {
      motif = m[1]
    } else {
      motif = info                      # secours : colonne déjà un motif brut
    }

    if (match(info, /ID=([^;]+)/, idm)) {
      name = idm[1]
    } else {
      name = info
    }

    print $1, start, $3, motif, name
  }
' "${RAW_BED}" > "${CLEAN_BED}"

echo "OK : ${OUTDIR}/${CLEAN_BED} ($(wc -l < "${CLEAN_BED}") loci)"
echo "Renseigne ce chemin dans config.yaml sous longtr.regions_bed"
