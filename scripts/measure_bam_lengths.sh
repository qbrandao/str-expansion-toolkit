#!/usr/bin/env bash
# Exact read count, total bases, depth, read length N50, ultra-long content and
# the spanning depth profile of one sample, from one or more BAM files.
#
# The spanning depth is the quantity that decides what a sample can measure. A
# read genotypes a repeat only if it covers the motif plus an anchor on each
# side, so for a locus of length W and an anchor f a read of length L can be
# placed spanning it in L - W - 2f positions, and the expected number of
# spanning reads is the sum of that over all reads divided by the genome
# length. Read length N50 on its own does not give this, and two samples at the
# same total depth can differ by orders of magnitude in it.
#
# All of it comes from one pass over the data, because a second pass over a
# cohort of raw BAMs is hours of I/O that buys nothing.
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
# Anchor required on each side of the repeat for a read to be counted as
# spanning it. 1000 bp is conservative and the value belongs in the methods.
ANCHOR="${ANCHOR:-1000}"
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

awk -v label="$LABEL" -v g="$GENOME_BP" -v nb="$#" -v anchor="$ANCHOR" '
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
    # spanning depth at the locus sizes that matter for repeat genotyping
    nw = split("1000 2000 5000 10000 20000 50000", W, " ")
    out = ""
    for (q = 1; q <= nw; q++) {
      w = W[q] + 0; sp = 0
      for (i = 1; i <= rows; i++) {
        eff = L[i] - w - 2 * anchor
        if (eff > 0) sp += eff * C[i]
      }
      out = out sprintf("\t%.3f", sp / g)
    }
    printf "%s\t%d\t%d\t%d\t%.3f\t%d\t%d\t%.0f\t%.3f\t%.3f\t%.4f\t%.4f\t%d%s\n",
      label, nb, n, tot, tot / g, n50, mx, tot / n,
      b100 / g, b50 / g, n100 / n, n50c / n, anchor, out
  }' "$WORK/hist.sorted" > "$WORK/row"

HDR="label	n_bam	n_reads	total_bases	depth_x	read_n50	max_len	mean_len	depth_over_100kb	depth_over_50kb	frac_reads_over_100kb	frac_reads_over_50kb	anchor_bp	span_1kb	span_2kb	span_5kb	span_10kb	span_20kb	span_50kb"
{ printf '%s\n' "$HDR"; cat "$WORK/row"; } > "$OUT.tmp.$$"
mv -f "$OUT.tmp.$$" "$OUT"
echo "written: $OUT"
printf '%s\n' "$HDR"; cat "$WORK/row"
