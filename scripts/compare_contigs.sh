#!/usr/bin/env bash
# Compare the contig set of one or more BAM files against a reference index, and
# against each other.
#
# Two BAMs are only comparable at a repeat locus if the chromosome holding it has
# the same name AND the same length in both. Equal lengths on the primary
# chromosomes mean the coordinates agree, so a catalogue built on one reference
# addresses the same bases in the other. A different accessory contig set, alts
# or decoys present on one side and absent on the other, does not move
# coordinates but does move ambiguous reads, which at repeat loci is a case
# control confounder rather than a cosmetic difference.
#
# Usage: compare_contigs.sh REFERENCE.fa[.fai] BAM [BAM...]
set -euo pipefail
REF="${1:?usage: compare_contigs.sh REFERENCE.fa[.fai] BAM [BAM...]}"
shift
command -v samtools >/dev/null || { echo "samtools not on PATH" >&2; exit 1; }

FAI="$REF"
[ "${FAI##*.}" = "fai" ] || FAI="$REF.fai"
if [ ! -s "$FAI" ]; then
  echo "no index at $FAI, build it with: samtools faidx $REF" >&2; exit 1
fi

W=$(mktemp -d "${TMPDIR:-/tmp}/ctg.XXXXXX")
trap 'rm -rf "$W"' EXIT

cut -f1,2 "$FAI" | LC_ALL=C sort > "$W/ref"
echo "reference $(basename "$REF"): $(wc -l < "$W/ref") contigs"
echo

PRIMARY='^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$'

i=0
for b in "$@"; do
  i=$((i + 1))
  [ -s "$b" ] || { echo "missing: $b" >&2; exit 1; }
  samtools view -H "$b" \
  | awk -F'\t' '$1=="@SQ"{ n=""; l=""
      for (j=2; j<=NF; j++) { if ($j ~ /^SN:/) n=substr($j,4); if ($j ~ /^LN:/) l=substr($j,4) }
      if (n != "") print n "\t" l }' | LC_ALL=C sort > "$W/b$i"
  tot=$(wc -l < "$W/b$i")
  same=$(LC_ALL=C comm -12 "$W/ref" "$W/b$i" | wc -l)
  onlyb=$(LC_ALL=C comm -13 "$W/ref" "$W/b$i" | wc -l)
  onlyr=$(LC_ALL=C comm -23 "$W/ref" "$W/b$i" | wc -l)
  # a name present in both but with a different length is the dangerous case
  clash=$(LC_ALL=C join -t'	' <(cut -f1,2 "$W/ref") <(cut -f1,2 "$W/b$i") \
          | awk -F'\t' '$2 != $3 {print}' | wc -l)
  pri_clash=$(LC_ALL=C join -t'	' <(cut -f1,2 "$W/ref") <(cut -f1,2 "$W/b$i") \
          | awk -F'\t' -v p="$PRIMARY" '$2 != $3 && $1 ~ p {print}' | wc -l)

  printf '%s\n' "$(basename "$b")"
  printf '   contigs in the BAM            %6d\n' "$tot"
  printf '   identical name and length     %6d\n' "$same"
  printf '   in the BAM, not the reference  %6d\n' "$onlyb"
  printf '   in the reference, not the BAM  %6d\n' "$onlyr"
  printf '   same name, DIFFERENT length    %6d   (of which primary chromosomes: %d)\n' "$clash" "$pri_clash"
  if [ "$pri_clash" -gt 0 ]; then
    echo "   VERDICT: incompatible, a primary chromosome differs in length, so the"
    echo "            coordinates do not agree and no catalogue is transferable."
    LC_ALL=C join -t'	' <(cut -f1,2 "$W/ref") <(cut -f1,2 "$W/b$i") \
      | awk -F'\t' -v p="$PRIMARY" '$2 != $3 && $1 ~ p {printf "            %s ref=%s bam=%s\n", $1, $2, $3}'
  elif [ "$onlyb" -gt 0 ] || [ "$onlyr" -gt 0 ]; then
    echo "   VERDICT: coordinates agree on the primary chromosomes, but the accessory"
    echo "            contig set differs. Mapping of ambiguous reads differs, so this"
    echo "            must be identical across cases and controls."
  else
    echo "   VERDICT: identical contig set."
  fi
  echo
done

if [ "$#" -gt 1 ]; then
  echo "=== the BAMs against each other ==="
  j=1
  while [ "$j" -lt "$#" ]; do
    k=$((j + 1))
    while [ "$k" -le "$#" ]; do
      if cmp -s "$W/b$j" "$W/b$k"; then
        printf '   BAM %d and %d: identical contig set\n' "$j" "$k"
      else
        printf '   BAM %d and %d: DIFFER (%d contigs only in %d, %d only in %d)\n' \
          "$j" "$k" \
          "$(LC_ALL=C comm -23 "$W/b$j" "$W/b$k" | wc -l)" "$j" \
          "$(LC_ALL=C comm -13 "$W/b$j" "$W/b$k" | wc -l)" "$k"
      fi
      k=$((k + 1))
    done
    j=$((j + 1))
  done
  echo
  echo "   (BAM numbering follows the order given on the command line)"
fi
