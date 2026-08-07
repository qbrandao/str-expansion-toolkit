# SLURM test scripts

Run from the repo root (`sbatch scripts/00_unit_tests.sbatch`, etc.), or
chain everything at once with `./scripts/submit_full_pipeline.sh`.

## Prerequisites (place at the repo root before running)

- `config.yaml` (copied from `config.example.yaml`, paths adapted)
- `controls.tsv`, `patients.tsv`: TSV with header `sample_id\tbam_path\tfastq_path`
- `genes.bed.gz`, `MANE_Select_exons.bed.gz`: used by `compare`
- `longtr.regions_bed` (config.yaml): see `./scripts/prepare_longtr_bed.sh`
  below to generate it from the Adotto catalog (requires `gawk`)
- Adapt `<your_partition>` in each `.sbatch` file (your cluster's SLURM partition)
- `mkdir -p logs` (done automatically by the scripts, but the directory
  must be writable from compute nodes)

## Recommended order (first run, step by step)

1. `00_unit_tests.sbatch` — validates the package installation (no bioinformatics tool is run)
2. `01_smoke_test_single_tool.sbatch` — validates `config.yaml`/paths with LongTR only (fast)
3. `02_detect_single_patient.sbatch` — validates the full pipeline (default tools) on 1 patient
4. `03_detect_controls_array.sbatch` — detection across the whole control cohort (job array,
   adjust `--array=1-N%4` where N = number of lines in `controls.tsv`)
5. `04_detect_patients_array.sbatch` — detection across all patients (job array, same idea)
6. `05_build_controls.sbatch` — after step 4 has fully completed (control JSON registry)
7. `06_compare.sbatch` — after steps 4 and 5 have completed (final report)

## Preparing the LongTR regions catalog

```bash
./scripts/prepare_longtr_bed.sh /mnt/references
```

Downloads the Adotto catalog (hg38) and converts it to the format expected
by LongTR (`chrom, start 1-based, end, motif[,motif2], name`). Requires
`gawk` (the script checks for it and exits otherwise).

## Chaining everything automatically

```bash
./scripts/submit_full_pipeline.sh
```

Submits 00 → 01 → (03 array + 04 array in parallel) → 05 → 06, with
`--dependency=afterok:...` computed automatically (job array size counted
dynamically from `controls.tsv`/`patients.tsv`). Track progress with
`squeue -u $USER`.
