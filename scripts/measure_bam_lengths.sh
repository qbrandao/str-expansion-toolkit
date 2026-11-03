#!/usr/bin/env bash
# Exact read count, total bases, depth, read length N50 and ultra-long content
# of one sample, from one or more BAM files.
#
# Primary alignments only (-F 0x900), so secondary and supplementary records
# cannot count a base twice. Works on unaligned run BAMs and on aligned BAMs.
#
# The N50 is computed from a length histogram sorted by an external sort rather
# than inside awk, because array sorting is a gawk extension and this has to
# run under mawk and busybox awk as well.
#
# Usage: measure_bam_lengths.sh OUT.tsv LABEL BAM [BAM...]
set -euo pipefail

[ "$#" -ge 3 ] || { echo "usage: $0 OUT.tsv LABEL BAM [BAM...]" >&2; exit 2; }
OUT="$1"; LABEL="$2"; shift 2

# GRCh38 primary assembly, chr1-22 X Y, no alts and no decoys. Override with
# GENOME_BP if you count coverage against a different span.
GENOME_BP="${GENOME_BP:-3088269832}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-6}}"

command -v samtools >/dev/null || { echo "samtools not on PATH" >&2; exit 1; }
for b in "$@"; do [ -s "$b" ] || { echo "missing or empty: $b" >&2; exit 1; }; done

WORK=$(mktemp -d "${TMPDIR:-/tmp}/measlen.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

# Pass 1: one line per distinct read length, "length count".
{
  set -o pipefail
  for b in "$@"; do
    samtools view -@ "$THREADS" -F 0x900 "$b" | awk '{print length($10)}'
  done
} | awk '{h[$1]++} END{for (l in h) print l, h[l]}' > "$WORK/hist"

[ -s "$WORK/hist" ] || { echo "no primary reads found in the given BAM" >&2; exit 1; }

# Pass 2: descending by length, accumulate bases until half the total.
sort -k1,1nr "$WORK/hist" > "$WORK/hist.sorted"

awk -v label="$LABEL" -v g="$GENOME_BP" -v nb="$#" '
  { len = $1 + 0; cnt = $2 + 0
    n += cnt; tot += len * cnt
    if (len > mx) mx = len
    if (len >= 100000) { n100 += cnt; b100 += len * cnt }
    if (len >=  50000) { n50c += cnt; b50  += len * cnt }
    L[NR] = len; C[NR] = cnt; rows = NR }
  END {
    if (n == 0) exit 1
    half = tot / 2; acc = 0; n50 = 0
    for (i = 1; i <= rows; i++) {
      acc += L[i] * C[i]
      if (acc >= half) { n50 = L[i]; break }
    }
    printf "%s\t%d\t%d\t%d\t%.3f\t%d\t%d\t%.0f\t%.3f\t%.3f\t%.4f\t%.4f\n",
      label, nb, n, tot, tot / g, n50, mx, tot / n,
      b100 / g, b50 / g, n100 / n, n50c / n
  }' "$WORK/hist.sorted" > "$WORK/row"

HDR="label	n_bam	n_reads	total_bases	depth_x	read_n50	max_len	mean_len	depth_over_100kb	depth_over_50kb	frac_reads_over_100kb	frac_reads_over_50kb"
{ printf '%s\n' "$HDR"; cat "$WORK/row"; } > "$OUT.tmp.$$"
mv -f "$OUT.tmp.$$" "$OUT"
echo "written: $OUT"
printf '%s\n' "$HDR"; cat "$WORK/row"
