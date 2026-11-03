#!/bin/bash
#SBATCH --job-name=str-merge-align
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --partition=<your_partition>
#SBATCH --array=1-N%8
# ^ Replace N with the number of lines in controls_samples.txt. The %8 throttle
#   is a disk constraint rather than a CPU one: each concurrent task holds its
#   own sorted output plus samtools sort temporaries.

# Merge the per-run unaligned BAMs of one control, subsample it to the common
# target depth and align it, producing one coordinate sorted, indexed BAM.
#
# The per-run reads are streamed through samtools fastq into minimap2, so no
# intermediate FASTQ is ever written. Methylation tags are carried as FASTQ
# comments by samtools fastq -T and restored by minimap2 -y.
#
# Subsampling happens before alignment, which is the whole point: the discarded
# reads are never aligned, and the LAST stage downstream, which dominates the
# runtime of the pipeline, sees only the retained reads. Random subsampling
# preserves the read length distribution, so the long reads that span an
# expansion are retained in proportion.
#
# Safety properties, matching the rest of the pipeline:
#   - pipefail, so a stage killed mid pipeline fails the step
#   - the output is renamed from a temporary name only after a clean exit, so a
#     truncated BAM is never accepted as finished
#   - an existing indexed output is left untouched, which makes the array
#     restartable after a preemption
#
# Prerequisites at the repo root: controls_samples.txt, downsample_plan.tsv
# (from scripts/make_downsample_plan.py) and the minimap2 index.

set -euo pipefail
mkdir -p logs

cd "${SLURM_SUBMIT_DIR}"

CONTROLS_ROOT="${CONTROLS_ROOT:-/store/brandaoq/controls}"
SAMPLES="controls_samples.txt"
PLAN="downsample_plan.tsv"
REF_MMI="${REF_MMI:-/store/brandaoq/ref/GRCh38.mmi}"
SUBSAMPLE_SEED="${SUBSAMPLE_SEED:-42}"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
SORT_THREADS=$(( THREADS / 4 )); [ "${SORT_THREADS}" -lt 1 ] && SORT_THREADS=1
SORT_MEM="2G"

SID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLES}")
[ -n "${SID}" ] || { echo "no sample at line ${SLURM_ARRAY_TASK_ID} of ${SAMPLES}" >&2; exit 1; }

DIR="${CONTROLS_ROOT}/${SID}"
OUT="${DIR}/${SID}.sorted.bam"
TMP="${DIR}/.${SID}.sorted.bam.tmp.$$"
TMPDIR_SORT="${SLURM_TMPDIR:-${DIR}}/sort_${SID}.$$"

echo "Control ${SLURM_ARRAY_TASK_ID}: ${SID} on $(hostname), ${THREADS} threads"

if [ -s "${OUT}" ] && { [ -s "${OUT}.bai" ] || [ -s "${OUT}.csi" ]; }; then
  echo "${SID} already aligned and indexed, nothing to do"; exit 0
fi

mapfile -t RUNS < <(find "${DIR}" -maxdepth 1 -type f -name '*.bam' \
                      ! -name "${SID}.sorted.bam" ! -name '.*' -printf '%p\n' | sort)
[ "${#RUNS[@]}" -gt 0 ] || { echo "no run BAM in ${DIR}" >&2; exit 1; }
[ -s "${REF_MMI}" ] || { echo "reference index not found: ${REF_MMI}" >&2; exit 1; }
[ -s "${PLAN}" ]    || { echo "downsampling plan not found: ${PLAN}" >&2; exit 1; }

# The fraction column is located by its header name rather than by position, so
# a change in the plan's layout cannot silently feed the wrong number in. A
# missing row is a hard error rather than a fall back to the native depth,
# because a single control left at full depth would be the deepest in the panel
# and would bias every depth dependent statistic computed from it.
FRAC=$(awk -F'\t' -v s="${SID}" '
  NR==1 { for (i=1; i<=NF; i++) if ($i=="subsample_fraction") col=i
          if (!col) { print "plan has no subsample_fraction column" > "/dev/stderr"; exit 2 }
          next }
  $1==s { print $col; found=1; exit }
  END   { if (!found) exit 1 }' "${PLAN}") \
  || { echo "no usable row for ${SID} in ${PLAN}" >&2; exit 1; }
case "${FRAC}" in
  ''|*[!0-9.]*) echo "unusable fraction for ${SID}: '${FRAC}'" >&2; exit 1 ;;
esac
awk -v f="${FRAC}" 'BEGIN{exit !(f>0 && f<=1)}' \
  || { echo "fraction out of range for ${SID}: ${FRAC}" >&2; exit 1; }
echo "${#RUNS[@]} run BAM to merge, subsampling fraction ${FRAC}"

mkdir -p "${TMPDIR_SORT}"
cleanup() { rm -rf "${TMPDIR_SORT}" "${TMP}"; }
trap cleanup EXIT

# One reader per run BAM, concatenated into a single minimap2 invocation, so the
# run boundaries disappear and the aligner sees one stream of reads.
{
  set -o pipefail
  for b in "${RUNS[@]}"; do
    if awk -v f="${FRAC}" 'BEGIN{exit !(f>=1)}'; then
      samtools fastq -@ 2 -T MM,ML "${b}"
    else
      samtools view -@ 2 -u --subsample "${FRAC}" --subsample-seed "${SUBSAMPLE_SEED}" "${b}" \
        | samtools fastq -@ 2 -T MM,ML -
    fi
  done
} | minimap2 -ax map-ont -Y -y --MD -t "${THREADS}" "${REF_MMI}" - \
  | samtools sort -@ "${SORT_THREADS}" -m "${SORT_MEM}" -T "${TMPDIR_SORT}/s" -o "${TMP}" -

samtools index -@ "${THREADS}" "${TMP}"
mv -f "${TMP}.bai" "${OUT}.bai"
mv -f "${TMP}" "${OUT}"
trap - EXIT; rm -rf "${TMPDIR_SORT}"

echo "${SID} done, $(stat -c %s "${OUT}") bytes"
samtools quickcheck -v "${OUT}" && echo "quickcheck ok"
