#!/usr/bin/env bash
# Spanning depth as a function of repeat length, for one sample.
#
# A read spans a locus only if it covers the repeat plus an anchor on each
# side, so the read length N50 on its own does not say what is callable. For a
# locus of length W and an anchor f, a read of length L can be placed so that
# it spans the locus in (L - W - 2f) positions, so the depth of spanning reads
# is sum over reads of max(0, L - W - 2f) divided by the genome length. That is
# the quantity that decides the largest expansion a sample can measure, and it
# is what this prints, per locus size.
#
# Run it on a patient and on a control and compare the two tables: the size at
# which the patient's spanning depth falls below the minimum read support used
# by the somatic analysis is the ceiling of the whole study, and loci above it
# are not comparable between the two.
#
# Usage: spanning_profile.sh BAM [BAM...]
#   ANCHOR=1000      anchor required on each side, in bp
#   GENOME_BP        defaults to GRCh38 primary, chr1-22 X Y
#   THREADS
set -euo pipefail

[ "$#" -ge 1 ] || { echo "usage: $0 BAM [BAM...]" >&2; exit 2; }
ANCHOR="${ANCHOR:-1000}"
GENOME_BP="${GENOME_BP:-3088269832}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-6}}"

command -v samtools >/dev/null || { echo "samtools not on PATH" >&2; exit 1; }
for b in "$@"; do [ -s "$b" ] || { echo "missing or empty: $b" >&2; exit 1; }; done

WORK=$(mktemp -d "${TMPDIR:-/tmp}/spanprof.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

{
  set -o pipefail
  for b in "$@"; do
    samtools view -@ "$THREADS" -F 0x900 "$b" | awk '{print length($10)}'
  done
} | awk '{h[$1]++} END{for (l in h) print l, h[l]}' > "$WORK/hist"

[ -s "$WORK/hist" ] || { echo "no primary reads found" >&2; exit 1; }

awk -v anchor="$ANCHOR" -v g="$GENOME_BP" '
  { len = $1 + 0; cnt = $2 + 0; L[++k] = len; C[k] = cnt; n += cnt; tot += len * cnt }
  END {
    split("0 100 200 500 1000 2000 3000 5000 8000 10000 15000 20000 30000 50000 100000", W, " ")
    printf "reads %d   bases %d   depth %.2fx   anchor %d bp each side\n\n", n, tot, tot/g, anchor
    printf "%12s  %14s  %10s\n", "locus_bp", "spanning_depth", "pct_of_raw"
    for (i = 1; i in W; i++) {
      w = W[i] + 0; s = 0
      for (j = 1; j <= k; j++) {
        eff = L[j] - w - 2 * anchor
        if (eff > 0) s += eff * C[j]
      }
      printf "%12d  %13.2fx  %9.1f%%\n", w, s / g, (tot ? 100 * s / tot : 0)
    }
  }' "$WORK/hist"
