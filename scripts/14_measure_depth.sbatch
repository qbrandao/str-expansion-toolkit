#!/bin/bash
#SBATCH --job-name=str-measure-depth
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --partition=<your_partition>
#SBATCH --array=1-N
# ^ Replace N with the number of lines in controls_samples.txt:
#     wc -l < controls_samples.txt
#   No concurrency throttle is needed: each task is I/O bound on a different
#   file and uses little CPU and memory.

# Measure the exact depth and read length N50 of every control, one array task
# per sample, by streaming each run BAM once.
#
# This exists because depth cannot be inferred from file size. The published
# HPRC release 2 index describes the files of that release rather than the files
# actually downloaded, so a calibration fitted against it has no predictive
# value. Measuring is one pass over the raw data, a few hundred CPU hours for a
# whole panel, which is negligible against genotyping and removes every
# estimated depth from the study.
#
# The N50 obtained here is also what decides whether the controls and the
# patients are comparable in read length, which subsampling cannot fix because
# removing reads at random preserves the length distribution.
#
# Each task writes its own small TSV, so no two tasks write the same file.
# Collect them once the array has finished:
#   { head -1 measured/*.tsv | head -1; tail -q -n +2 measured/*.tsv; } \
#     | sort -u > measured_depth.tsv

set -euo pipefail
mkdir -p logs

cd "${SLURM_SUBMIT_DIR}"

CONTROLS_ROOT="${CONTROLS_ROOT:-/store/brandaoq/controls}"
SAMPLES="controls_samples.txt"
OUTDIR="measured"

mkdir -p "${OUTDIR}"
SID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLES}")
[ -n "${SID}" ] || { echo "no sample at line ${SLURM_ARRAY_TASK_ID} of ${SAMPLES}" >&2; exit 1; }

OUT="${OUTDIR}/${SID}.tsv"
if [ -s "${OUT}" ]; then
  echo "${SID} already measured, skipping"; exit 0
fi

DIR="${CONTROLS_ROOT}/${SID}"
mapfile -t RUNS < <(find "${DIR}" -maxdepth 1 -type f -name '*.bam' \
                      ! -name "${SID}.sorted.bam" ! -name '.*' -printf '%p\n' | sort)
[ "${#RUNS[@]}" -gt 0 ] || { echo "no run BAM in ${DIR}" >&2; exit 1; }

# A sample whose runs do not share one basecalling model must not be measured as
# a single unit, because it will not be merged as a single unit either. Resolve
# it first, by removing the odd run or by excluding the sample.
MODELS=$(printf '%s\n' "${RUNS[@]##*/}" \
         | grep -oE '_(sup|hac|fast)[0-9]+\.[0-9]+\.[0-9]+' | sort -u | paste -sd, -)
case "${MODELS}" in
  *,*) echo "${SID} mixes basecalling models (${MODELS}), resolve it before measuring" >&2
       exit 1 ;;
esac

echo "Control ${SLURM_ARRAY_TASK_ID}: ${SID}, ${#RUNS[@]} run BAM, model ${MODELS:-unknown}"
THREADS="${SLURM_CPUS_PER_TASK}" ./scripts/measure_bam_lengths.sh "${OUT}" "${SID}" "${RUNS[@]}"
