#!/bin/bash
# Submits the full pipeline with the correct SLURM dependencies:
#
#   00 unit tests
#        |
#   01 smoke test (1 tool)
#        |
#        +--> 03 detect controls (array) --+--> 05 build-controls --> 06 compare
#        |                                 |                             ^
#        |                                 +--> 07 repertoire            |
#        +--> 04 detect patients (array) ------------------------------- +
#        |
#        +--> 08 detect CEPH (array) --> 09 meiotic-instability
#
#   10 somatic-instability is NOT chained here: it needs all sample
#   directories under a single --detect-dir, which depends on how you lay
#   out controls/patients/CEPH outputs. Submit it manually once the relevant
#   detect jobs are done (see the note in 10_somatic_instability.sbatch).
#
# Usage: ./scripts/submit_full_pipeline.sh
# Run from the repo root (where config.yaml lives).

set -euo pipefail
mkdir -p logs

n_lines() { tail -n +2 "$1" | wc -l; }

N_CONTROLS=$(n_lines controls.tsv)
N_PATIENTS=$(n_lines patients.tsv)
N_CEPH=$(n_lines ceph.tsv)

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

echo "== 08: detect CEPH quartet (array, ${N_CEPH} samples) =="
JID_08=$(sbatch --parsable --dependency=afterok:${JID_01} --array=1-${N_CEPH}%2 scripts/08_detect_ceph_array.sbatch)
echo "  job ${JID_08}"

echo "== 05: build-controls (after the full 03 array) =="
JID_05=$(sbatch --parsable --dependency=afterok:${JID_03} scripts/05_build_controls.sbatch)
echo "  job ${JID_05}"

echo "== 07: repertoire (after the full 03 array) =="
JID_07=$(sbatch --parsable --dependency=afterok:${JID_03} scripts/07_repertoire.sbatch)
echo "  job ${JID_07}"

echo "== 06: compare (after the 04 array and job 05) =="
JID_06=$(sbatch --parsable --dependency=afterok:${JID_04},${JID_05} scripts/06_compare.sbatch)
echo "  job ${JID_06}"

echo "== 09: meiotic-instability (after the full 08 array) =="
JID_09=$(sbatch --parsable --dependency=afterok:${JID_08} scripts/09_meiotic_instability.sbatch)
echo "  job ${JID_09}"

echo
echo "Pipeline submitted. Track progress: squeue -u \$USER"
echo
echo "Expected outputs:"
echo "  results/report.tsv                      (job ${JID_06}) -- patient vs control comparison"
echo "  results/repertoire.tsv                  (job ${JID_07}) -- genome-wide VNTR repertoire"
echo "  results/repertoire_summary.tsv          (job ${JID_07}) -- counts by location x motif"
echo "  results/meiotic_instability.tsv         (job ${JID_09}) -- per-locus germline instability"
echo "  results/meiotic_instability_by_duo.tsv  (job ${JID_09}) -- per-duo summary"
echo
echo "Then submit somatic instability manually once sample directories are laid out:"
echo "  sbatch scripts/10_somatic_instability.sbatch"
