#!/bin/bash
# Triage of an input BAM before running the pipeline on it.
#
# Answers the questions that decide whether this BAM can be used as is, and
# with which configuration:
#   1. is it a long-read BAM at all, and is it readable and coordinate-sorted
#   2. which reference build and contig naming, and does it match config.yaml
#   3. which basecaller and chemistry (R9.4.1 vs R10.4.1), which drives the
#      clair3 model and the minimap2 preset
#   4. observed error rate, as a chemistry check independent of the header
#   5. depth, which drives LongTR --min-reads
#   6. whole genome or targeted, which decides whether a genome-wide
#      repertoire run makes sense
#   7. read groups and sample names
#   8. methylation tags
#
# Deliberately cheap: it reads the header, the index summary, and a bounded
# sample of reads from a few small regions. It is safe to run on a login
# node and takes well under a minute. It writes nothing except an optional
# report file.
#
# Usage:
#   ./inspect_bam.sh /path/to/sample.bam [/path/to/reference.fa] > report.txt

set -uo pipefail

BAM="${1:-}"
REF="${2:-}"

if [[ -z "${BAM}" ]]; then
  echo "usage: $0 <sample.bam> [reference.fa]" >&2
  exit 2
fi
if [[ ! -f "${BAM}" ]]; then
  echo "ERROR: not found: ${BAM}" >&2
  exit 2
fi
command -v samtools >/dev/null || { echo "ERROR: samtools not on PATH" >&2; exit 2; }

hr() { printf '%s\n' "----------------------------------------------------------------"; }
sec() { hr; printf '%s\n' "$1"; hr; }

echo "BAM      : ${BAM}"
echo "size     : $(du -h "${BAM}" | cut -f1)"
echo "samtools : $(samtools --version | head -1)"
echo "date     : $(date -Is)"
echo

# --- 1. integrity, sort order, index ---------------------------------
sec "1. INTEGRITY, SORT ORDER, INDEX"
if samtools quickcheck -v "${BAM}"; then
  echo "quickcheck: OK (complete file, EOF block present)"
else
  echo "quickcheck: FAILED. The file is truncated or corrupt. Do not use it."
fi

SO="$(samtools view -H "${BAM}" | awk '/^@HD/{for(i=1;i<=NF;i++) if($i ~ /^SO:/) print substr($i,4)}')"
echo "sort order (@HD SO): ${SO:-absent}"
if [[ "${SO}" != "coordinate" ]]; then
  echo "  NOTE: every tool here needs a coordinate-sorted BAM."
  echo "        Sort it first: samtools sort -@ 8 -o sorted.bam ${BAM}"
fi

IDX="none"
for c in "${BAM}.bai" "${BAM}.csi" "${BAM%.bam}.bai"; do
  [[ -f "${c}" ]] && IDX="${c}"
done
echo "index: ${IDX}"
if [[ "${IDX}" == "none" ]]; then
  echo "  NOTE: no index. str-toolkit will symlink the BAM into the output"
  echo "        directory and index the symlink, leaving your directory untouched."
fi

# --- 2. reference build and contig naming ----------------------------
sec "2. REFERENCE BUILD AND CONTIG NAMING"
samtools view -H "${BAM}" | awk '/^@SQ/{n++} END{print "contigs in header: " n+0}'
echo "first 5 contigs:"
samtools view -H "${BAM}" | awk '/^@SQ/{for(i=1;i<=NF;i++){if($i ~ /^SN:/) sn=substr($i,4); if($i ~ /^LN:/) ln=substr($i,4)} print "  " sn "\t" ln}' | head -5

CHR1LEN="$(samtools view -H "${BAM}" | awk '/^@SQ/{sn="";ln="";for(i=1;i<=NF;i++){if($i ~ /^SN:/) sn=substr($i,4); if($i ~ /^LN:/) ln=substr($i,4)} if(sn=="chr1"||sn=="1") print ln}' | head -1)"
case "${CHR1LEN}" in
  248956422) echo "build: GRCh38 / hg38 (chr1 length ${CHR1LEN})" ;;
  249250621) echo "build: GRCh37 / hg19 (chr1 length ${CHR1LEN})" ;;
  248387328) echo "build: T2T-CHM13v2.0 (chr1 length ${CHR1LEN})" ;;
  "")        echo "build: could not identify (no chr1 or 1 contig in the header)" ;;
  *)         echo "build: unrecognised, chr1 length ${CHR1LEN}" ;;
esac

NAMING="$(samtools view -H "${BAM}" | awk '/^@SQ/{for(i=1;i<=NF;i++) if($i ~ /^SN:/) print substr($i,4)}' | grep -qx 'chr1' && echo "chr-prefixed (chr1, chr2, ...)" || echo "unprefixed (1, 2, ...)")"
echo "contig naming: ${NAMING}"
echo "  NOTE: the catalogs (VAMOS, Adotto/LongTR, TRF) must use the same"
echo "        naming, otherwise every locus silently fails to match."

if [[ -n "${REF}" && -f "${REF}" ]]; then
  echo
  echo "comparing against the reference given on the command line:"
  if [[ -f "${REF}.fai" ]]; then
    samtools view -H "${BAM}" | awk '/^@SQ/{sn="";ln="";for(i=1;i<=NF;i++){if($i ~ /^SN:/) sn=substr($i,4); if($i ~ /^LN:/) ln=substr($i,4)} print sn"\t"ln}' | sort > /tmp/.bam_sq.$$
    awk '{print $1"\t"$2}' "${REF}.fai" | sort > /tmp/.ref_sq.$$
    NDIFF="$(comm -3 /tmp/.bam_sq.$$ /tmp/.ref_sq.$$ | wc -l)"
    if [[ "${NDIFF}" -eq 0 ]]; then
      echo "  contigs and lengths are identical. The BAM can be reused as is."
    else
      echo "  ${NDIFF} contig line(s) differ. Showing up to 10:"
      comm -3 /tmp/.bam_sq.$$ /tmp/.ref_sq.$$ | head -10 | sed 's/^/    /'
      echo "  If the shared contigs match in length, the difference is only in the"
      echo "  extra contigs and the BAM is still usable. If chr1..chr22 lengths"
      echo "  differ, this is another build: rerun detect with --realign."
    fi
    rm -f /tmp/.bam_sq.$$ /tmp/.ref_sq.$$
  else
    echo "  ${REF}.fai is missing, run: samtools faidx ${REF}"
  fi
fi

# --- 3. basecaller and chemistry from the header ---------------------
sec "3. BASECALLER AND CHEMISTRY (from the header)"
echo "@RG lines:"
samtools view -H "${BAM}" | grep '^@RG' | head -10 | sed 's/^/  /' || echo "  none"
echo
echo "@PG program chain:"
samtools view -H "${BAM}" | grep '^@PG' | sed -E 's/\t/ /g' | cut -c1-300 | head -12 | sed 's/^/  /' || echo "  none"
echo
HDR="$(samtools view -H "${BAM}")"
echo "chemistry markers found in the header:"
for pat in r9.4.1 r941 r10.4.1 r1041 e8.2 e82 dna_r9 dna_r10 guppy dorado bonito 400bps 260bps 450bps _sup _hac _fast; do
  n="$(grep -o -i -- "${pat}" <<<"${HDR}" | wc -l)"
  [[ "${n}" -gt 0 ]] && printf '  %-10s %s\n' "${pat}" "x${n}"
done
if ! grep -q -i -E 'r9|r10|guppy|dorado|bonito' <<<"${HDR}"; then
  echo "  none. The header carries no basecaller information, so rely on the"
  echo "  error rate measured in section 4 instead."
fi

# --- 4. observed error rate and read lengths -------------------------
sec "4. OBSERVED ERROR RATE AND READ LENGTHS (20000 primary reads)"
echo "Sampled from a few autosomal windows, primary alignments only."
echo "Reference: R9.4.1 sits around 3 to 6 percent, R10.4.1 sup around 0.5 to 2 percent."
echo
REGIONS=""
if grep -q 'SN:chr1' <<<"${HDR}"; then
  REGIONS="chr1:20000000-21000000 chr7:50000000-51000000 chr12:60000000-61000000"
elif grep -q 'SN:1' <<<"${HDR}"; then
  REGIONS="1:20000000-21000000 7:50000000-51000000 12:60000000-61000000"
fi

if [[ -z "${REGIONS}" || "${IDX}" == "none" ]]; then
  echo "(no index or unknown contig names: sampling the first reads of the file instead)"
  SAMPLE_CMD="samtools view -F 0x904 ${BAM}"
else
  SAMPLE_CMD="samtools view -F 0x904 ${BAM} ${REGIONS}"
fi

# One pass writes three numbers per read (read length, NM, alignment length
# from the CIGAR) to a temporary file. The N50 is then taken from a sort,
# rather than from awk's asort, which only exists in gawk.
TMPLEN="$(mktemp)"
trap 'rm -f "${TMPLEN}"' EXIT

${SAMPLE_CMD} 2>/dev/null | head -20000 | awk '
{
  nm = -1
  for (i = 12; i <= NF; i++) if ($i ~ /^NM:i:/) nm = substr($i, 6) + 0
  cig = $6
  alen = 0
  n = ""
  for (j = 1; j <= length(cig); j++) {
    c = substr(cig, j, 1)
    if (c ~ /[0-9]/) n = n c
    else { if (c ~ /[MIDX=]/) alen += n + 0; n = "" }
  }
  rl = length($10)
  if (rl > 0) print rl "\t" nm "\t" alen
}' > "${TMPLEN}"

NSAMPLED="$(wc -l < "${TMPLEN}")"
if [[ "${NSAMPLED}" -eq 0 ]]; then
  echo "  no primary read sampled. Is this BAM empty, or unaligned?"
else
  TOTLEN="$(awk -F'\t' '{s += $1} END {print s}' "${TMPLEN}")"
  N50="$(cut -f1 "${TMPLEN}" | sort -nr | awk -v tot="${TOTLEN}" \
        '{ acc += $1; if (acc >= tot / 2) { print $1; exit } }')"

  awk -F'\t' -v n50="${N50}" '
  { nreads++; sumlen += $1; if ($1 > maxlen) maxlen = $1
    if ($2 >= 0 && $3 > 0) { sumnm += $2; sumalen += $3; nerr++ } }
  END {
    mean = sumlen / nreads
    printf "  reads sampled          : %d\n", nreads
    printf "  mean read length       : %.0f bp\n", mean
    printf "  read length N50        : %d bp\n", n50
    printf "  max read length        : %d bp\n", maxlen
    if (nerr > 0) printf "  error rate (NM/alnlen) : %.2f %%\n", 100 * sumnm / sumalen
    else          print  "  error rate             : no NM tag on the sampled reads, cannot be measured"
    print ""
    # Read length settles whether this is long-read data at all, so it is
    # checked before the chemistry: a short-read BAM has a low error rate
    # too, and would otherwise be reported as good R10 data.
    if (mean < 500) {
      print "  --> mean read length under 500 bp. This is not long-read data, and"
      print "      none of these tools can be used on it. Stop here."
    } else if (nerr == 0) {
      print "  --> chemistry cannot be judged without NM tags. Rely on section 3."
    } else {
      err = 100 * sumnm / sumalen
      if (err >= 3.0)      print "  --> consistent with R9.4.1. Use an R9 clair3 model, not r1041_e82_400bps_sup_v500."
      else if (err >= 2.0) print "  --> between the two chemistries. Settle it from section 3 before choosing the clair3 model."
      else                 print "  --> consistent with R10.4.1 sup. r1041_e82_400bps_sup_v500 is the right clair3 model."
    }
  }' "${TMPLEN}"
fi

# --- 5. depth ---------------------------------------------------------
sec "5. DEPTH"
if [[ "${IDX}" != "none" && -n "${REGIONS}" ]]; then
  echo "samtools coverage on the sampled windows:"
  # shellcheck disable=SC2086
  for r in ${REGIONS}; do
    samtools coverage -r "${r}" "${BAM}" 2>/dev/null | tail -n +2 | \
      awk -v r="${r}" '{printf "  %-28s covered %6.2f %%   mean depth %6.2f x\n", r, $6, $7}'
  done
  echo
  echo "  NOTE: LongTR defaults to --min-reads 10. Below about 15x it will drop"
  echo "        many loci, and VAMOS phasing through clair3 becomes unreliable."
else
  echo "no index, so depth cannot be sampled by region."
fi

# --- 6. whole genome or targeted -------------------------------------
sec "6. WHOLE GENOME OR TARGETED"
if [[ "${IDX}" != "none" ]]; then
  samtools idxstats "${BAM}" 2>/dev/null | awk '
    $1 != "*" && $2 > 1000000 { n++; tot += $3; if ($3 > 0) withreads++ }
    END {
      printf "  main contigs (>1 Mb)        : %d\n", n
      printf "  of which carrying reads     : %d\n", withreads
      printf "  total mapped reads on those : %d\n", tot
      if (n > 0 && withreads < n * 0.5)
        print "  --> reads on fewer than half the contigs. This looks TARGETED, not whole genome."
      else
        print "  --> reads spread across the contigs, consistent with whole genome."
    }'
  echo
  echo "  top 5 contigs by mapped reads:"
  samtools idxstats "${BAM}" 2>/dev/null | awk '$1 != "*"' | sort -k3,3nr | head -5 | awk '{printf "    %-12s %12d reads\n", $1, $3}'
  echo
  echo "  unmapped reads (idxstats '*'): $(samtools idxstats "${BAM}" 2>/dev/null | awk '$1=="*"{print $4}')"
else
  echo "no index, so idxstats is unavailable."
fi

# --- 7. read groups and sample names ---------------------------------
sec "7. READ GROUPS AND SAMPLE NAMES"
SM="$(samtools view -H "${BAM}" | awk '/^@RG/{for(i=1;i<=NF;i++) if($i ~ /^SM:/) print substr($i,4)}' | sort -u)"
NRG="$(samtools view -H "${BAM}" | grep -c '^@RG' || true)"
echo "  @RG count   : ${NRG}"
echo "  SM value(s) : ${SM:-none}"
if [[ "$(wc -w <<<"${SM}")" -gt 1 ]]; then
  echo "  --> several samples in one BAM. Split it per sample before genotyping,"
  echo "      otherwise the genotypes are a mixture."
fi
if [[ "${NRG}" -eq 0 ]]; then
  echo "  --> no @RG. Harmless here: LongTR is called with --bam-samps/--bam-libs,"
  echo "      and whatshap with --ignore-read-groups."
fi

# --- 8. methylation tags ---------------------------------------------
sec "8. METHYLATION TAGS"
MMCOUNT="$(${SAMPLE_CMD} 2>/dev/null | head -2000 | grep -c 'MM:Z:' || true)"
echo "  reads carrying MM:Z among the first 2000 sampled: ${MMCOUNT}"
if [[ "${MMCOUNT}" -eq 0 ]]; then
  echo "  --> no methylation tags. Not needed for repeat sizing, but a methylation"
  echo "      analysis at expanded loci would not be possible from this BAM."
fi

# --- verdict ----------------------------------------------------------
sec "VERDICT AND NEXT COMMAND"
cat <<'EOT'
Read sections 1, 2 and 4 first. In short:

  section 1 not coordinate-sorted   -> sort before anything else
  section 2 build differs from your reference -> add --realign to detect,
            which needs a FASTQ, or realign the BAM yourself
  section 2 contig naming differs from the catalogs -> fix the catalogs,
            never the BAM
  section 4 error rate at or above 2.5 percent -> R9.4.1, so use an R9
            clair3 model rather than r1041_e82_400bps_sup_v500
  section 4 mean read length under 500 bp -> not long-read data, stop here
  section 6 targeted -> a genome-wide repertoire run is not meaningful,
            restrict the catalogs to the captured regions

If sections 1, 2, 4 and 6 are all fine, start with the fast single-tool test
(LongTR only, no clair3, no LAST), then the full three-tool run:

  sbatch scripts/12_test_new_bam_longtr.sbatch
  sbatch scripts/13_test_new_bam_full.sbatch
EOT
