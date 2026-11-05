#!/usr/bin/env bash
# Decide whether a read length N50 measured from primary alignments is the real
# read length of the library, or an artefact of how the BAM was built.
#
# Three checks, a header read plus two streaming passes.
#
# 1. The aligner command line from @PG. With minimap2 -Y every supplementary
#    alignment carries the full read sequence, which inflates the file without
#    changing the primary length distribution. Without -Y they are hard clipped.
# 2. A record census by flag class, with the ratio of all stored bases to
#    primary bases. A ratio well above 1 with -Y explains a file far larger than
#    its primary sequence content.
# 3. Whether every read name holds exactly one primary record. If so the
#    measured lengths are the real read lengths. If not, the reads were
#    fragmented and the distribution understates them.
#
# Flags are decomposed arithmetically rather than with and(), which is a gawk
# extension and not available under mawk.
#
# Usage: bam_length_sanity.sh BAM [THREADS]
set -euo pipefail
BAM="${1:?usage: bam_length_sanity.sh BAM [THREADS]}"
THREADS="${2:-6}"
[ -s "$BAM" ] || { echo "missing or empty: $BAM" >&2; exit 1; }
command -v samtools >/dev/null || { echo "samtools not on PATH" >&2; exit 1; }

echo "=== 1. aligner command line (@PG) ==="
samtools view -H "$BAM" | awk -F'\t' '$1=="@PG"{for(i=2;i<=NF;i++) print "   " $i}' || true
NY=$(samtools view -H "$BAM" | grep -c 'CL:.*minimap2.*-[a-zA-Z]*Y' || true)
echo "   @PG lines whose minimap2 command carries -Y : ${NY}"
echo

echo "=== 2. record census (one pass) ==="
samtools view -@ "$THREADS" "$BAM" \
| awk '{
    f = $2; L = length($10); n++; tb += L
    unm  = int(f /    4) % 2
    sec  = int(f /  256) % 2
    sup  = int(f / 2048) % 2
    if (unm)      { nu++ }
    else if (sec) { ns++ }
    else if (sup) { np2++ }
    else          { npri++; pb += L }
  }
  END {
    printf "   records total        %15d\n", n
    printf "   primary mapped       %15d\n", npri
    printf "   secondary            %15d\n", ns
    printf "   supplementary        %15d\n", np2
    printf "   unmapped             %15d\n", nu
    printf "   bases in primary     %15d\n", pb
    printf "   bases in all records %15d\n", tb
    if (pb > 0) printf "   stored / primary     %15.2f\n", tb / pb
    print npri > "/dev/stderr"
  }' 2> /tmp/.npri.$$
NPRI=$(cat /tmp/.npri.$$); rm -f /tmp/.npri.$$
echo

echo "=== 3. one primary record per read name? (one pass plus a sort) ==="
NNAME=$(samtools view -@ "$THREADS" -F 0x904 "$BAM" | cut -f1 \
        | LC_ALL=C sort -u -S 2G --parallel="$THREADS" | wc -l)
echo "   primary records      ${NPRI}"
echo "   distinct read names  ${NNAME}"
if [ "$NPRI" = "$NNAME" ]; then
  echo "   -> one primary per read name, the measured read lengths are real"
else
  echo "   -> $(( NPRI - NNAME )) extra primary records beyond the distinct names."
  echo "      The reads are fragmented and the length distribution understates them."
  echo "      Recompute it per read name before using it."
fi
