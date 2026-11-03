#!/usr/bin/env bash
# Inventory of the downloaded HPRC control BAMs.
# Reports, per sample, the number of run BAMs, the pore chemistry, the dorado
# version, the basecalling model and the total size, then flags anything
# inconsistent. Read only, writes a single TSV. Usage:
#   bash audit_controls.sh /store/brandaoq/controls [ids.txt]
set -u
ROOT="${1:-/store/brandaoq/controls}"
IDLIST="${2:-}"
OUT="$PWD/controls_disk_audit.tsv"

printf 'sample_id\tn_bam\tchemistry\tdorado\tmodel\tmods\ttotal_bytes\ttotal_gib\tflag\n' > "$OUT"

if [ -n "$IDLIST" ] && [ -f "$IDLIST" ]; then
  mapfile -t SAMPLES < <(grep -oE '\b(HG|NA)[0-9]{5}\b' "$IDLIST" | sort -u)
else
  mapfile -t SAMPLES < <(find "$ROOT" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)
fi

for s in "${SAMPLES[@]}"; do
  d="$ROOT/$s"
  if [ ! -d "$d" ]; then
    printf '%s\t0\t-\t-\t-\t-\t0\t0.00\tMISSING_DIR\n' "$s" >> "$OUT"; continue
  fi
  mapfile -t bams < <(find "$d" -maxdepth 1 -type f -name '*.bam' -printf '%f\n' | sort)
  n=${#bams[@]}
  if [ "$n" -eq 0 ]; then
    printf '%s\t0\t-\t-\t-\t-\t0\t0.00\tNO_BAM\n' "$s" >> "$OUT"; continue
  fi
  chem=$(printf '%s\n' "${bams[@]}" | grep -oE 'R[0-9]{3,4}' | sort -u | paste -sd, -)
  dor=$(printf '%s\n'  "${bams[@]}" | grep -oE '(dorado|guppy)[0-9]+\.[0-9]+\.[0-9]+' | sort -u | paste -sd, -)
  mdl=$(printf '%s\n'  "${bams[@]}" | grep -oE '_(sup|hac|fast)[0-9]+\.[0-9]+\.[0-9]+' | sed 's/^_//' | sort -u | paste -sd, -)
  mods=$(printf '%s\n' "${bams[@]}" | grep -oE '5mCG(_5hmCG)?' | sort -u | paste -sd, -)
  bytes=0
  for b in "${bams[@]}"; do bytes=$(( bytes + $(stat -c %s "$d/$b") )); done
  gib=$(awk -v b="$bytes" 'BEGIN{printf "%.2f", b/1073741824}')

  flag=ok
  [ -z "$chem" ] && flag=NO_CHEM_TOKEN
  case "$chem" in *,*) flag=MIXED_CHEMISTRY ;; esac
  case "$chem" in R941|*,R941|R941,*) flag=R9_PRESENT ;; esac
  case "$dor" in  *,*) flag="${flag}+MIXED_BASECALLER" ;; esac
  case "$dor" in  *guppy*) flag="${flag}+GUPPY_R9_ERA" ;; esac
  case "$mdl" in  *,*) flag="${flag}+MIXED_MODEL" ;; esac
  case "$mdl" in  *hac*|*fast*) flag="${flag}+NOT_SUP" ;; esac
  [ -z "$mods" ] && flag="${flag}+NO_METH_IN_NAME"

  printf '%s\t%d\t%s\t%s\t%s\t%s\t%d\t%s\t%s\n' \
    "$s" "$n" "${chem:--}" "${dor:--}" "${mdl:--}" "${mods:--}" "$bytes" "$gib" "$flag" >> "$OUT"
done

echo "written: $OUT"
echo
echo "== samples =="
awk -F'\t' 'NR>1{n++} END{print "  directories audited : " n}' "$OUT"
awk -F'\t' 'NR>1 && $9!~/MISSING_DIR|NO_BAM/{n++} END{print "  with at least one BAM: " n+0}' "$OUT"
echo
echo "== chemistry =="      ; awk -F'\t' 'NR>1{c[$3]++} END{for(k in c) printf "  %-12s %d\n", k, c[k]}' "$OUT"
echo "== basecaller ==" ; awk -F'\t' 'NR>1{c[$4]++} END{for(k in c) printf "  %-12s %d\n", k, c[k]}' "$OUT"
echo "== model =="          ; awk -F'\t' 'NR>1{c[$5]++} END{for(k in c) printf "  %-12s %d\n", k, c[k]}' "$OUT"
echo "== flags =="          ; awk -F'\t' 'NR>1{c[$9]++} END{for(k in c) printf "  %-28s %d\n", k, c[k]}' "$OUT"
echo
awk -F'\t' 'NR>1{b+=$7; n+=$2} END{printf "== totals ==\n  run BAMs: %d\n  raw uBAM: %.2f TiB\n", n, b/1099511627776}' "$OUT"
echo
echo "== BAMs per sample =="; awk -F'\t' 'NR>1 && $2>0{c[$2]++} END{for(k in c) printf "  %s BAM: %d samples\n", k, c[k]}' "$OUT" | sort -V
