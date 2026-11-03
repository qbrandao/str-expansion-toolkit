#!/usr/bin/env bash
# Print the run BAM names of the samples the audit flagged, so the three
# unparsed ones and the four heterogeneous ones can be resolved by eye.
set -u
TSV="${1:-controls_disk_audit.tsv}"
ROOT="${2:-/store/brandaoq/controls}"
awk -F'\t' 'NR>1 && ($9 ~ /NO_CHEM_TOKEN|MIXED_BASECALLER|MIXED_MODEL|MIXED_CHEMISTRY|NOT_SUP/){print $1"\t"$9}' "$TSV" \
| while IFS=$'\t' read -r s flag; do
    printf '\n== %s   [%s]\n' "$s" "$flag"
    find "$ROOT/$s" -maxdepth 1 -type f -name '*.bam' -printf '   %10s  %f\n' 2>/dev/null \
      | sort -k2 || echo "   (unreadable)"
  done
