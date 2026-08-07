#!/bin/bash
# Submits the full test pipeline with the correct SLURM dependencies:
#
#   00 unit tests
#   01 smoke test (1 tool)                 03 detect controls (array)
#                                                    |
#                                           05 build-controls (afterok:03)
#   04 detect patients (array)  ---------------------|
#                                                    |
#                                           06 compare (afterok:04,05)
#
# Usage: ./scripts/submit_full_pipeline.sh
# Run from the repo root (where config.yaml lives).

set -euo pipefail
mkdir -p logs

n_lines() { tail -n +2 "$1" | wc -l; }

N_CONTROLS=$(n_lines controls.tsv)
N_PATIENTS=$(n_lines patients.tsv)

echo "== 00: unit tests =="
JID_00=$(sbatch --parsable scripts/00_unit_tests.sbatch)
echo "  job ${JID_00}"

echo "== 01: smoke test (1 tool) =="
JID_01=$(sbatch --parsable --dependency=afterok:${JID_00} scripts/01_smoke_test_single_tool.sbatch)
echo "  job ${JID_01}"

echo "== 03: detect controls (array, ${N_CONTROLS} samples) =="
JID_03=$(sbatch --parsable --dependency=afterok:${JID_01} --array=1-${N_CONTROLS}%4 scripts/03_detect_controls_array.sbatch)
echo "  job ${JID_03}"

echo "== 04: detect patients (array, ${N_PATIENTS} samples) =="
JID_04=$(sbatch --parsable --dependency=afterok:${JID_01} --array=1-${N_PATIENTS}%4 scripts/04_detect_patients_array.sbatch)
echo "  job ${JID_04}"

echo "== 05: build-controls (after the full 03 array) =="
JID_05=$(sbatch --parsable --dependency=afterok:${JID_03} scripts/05_build_controls.sbatch)
echo "  job ${JID_05}"

echo "== 06: compare (after the 04 array and job 05) =="
JID_06=$(sbatch --parsable --dependency=afterok:${JID_04},${JID_05} scripts/06_compare.sbatch)
echo "  job ${JID_06}"

echo
echo "Pipeline submitted. Track progress: squeue -u \$USER"
echo "Final report expected in results/report.tsv once job ${JID_06} completes."
