# SLURM scripts

Run from the repo root (`sbatch scripts/00_unit_tests.sbatch`, etc.), or
chain everything at once with `./scripts/submit_full_pipeline.sh`.

## Prerequisites (place at the repo root before running)

- `config.yaml` (copied from `config.example.yaml`, paths adapted)
- `controls.tsv`, `patients.tsv`, `ceph.tsv`: TSV with header
  `sample_id\tbam_path\tfastq_path`
- `duos_ceph.tsv`: TSV with header `duo_id\tparent_id\tchild_id\tduo_type`
  (see "Germline instability" below)
- `all_samples.tsv`: TSV with at least a `sample_id` column, for somatic
  instability across the whole cohort
- `genes.bed.gz`, `MANE_Select_exons.bed.gz`: used by `compare`,
  `repertoire` and both instability analyses
- `longtr.regions_bed` (config.yaml): see `./scripts/prepare_longtr_bed.sh`
  below to generate it from the Adotto catalog (requires `gawk`)
- Adapt `<your_partition>` in each `.sbatch` file
- `mkdir -p logs` (done automatically by the scripts, but the directory
  must be writable from compute nodes)

## Scripts

| Script | Purpose | Resources |
|---|---|---|
| `00_unit_tests.sbatch` | pytest only, no bioinformatics tools | 2 CPU / 4G / 15min |
| `01_smoke_test_single_tool.sbatch` | `detect --tools longtr` on 1 sample, validates config | 8 CPU / 16G / 4h |
| `02_detect_single_patient.sbatch` | full `detect` on 1 patient | 16 CPU / 64G / 24h |
| `03_detect_controls_array.sbatch` | `detect` across the control cohort | idem, xN |
| `04_detect_patients_array.sbatch` | `detect` across all patients | idem, xN |
| `05_build_controls.sbatch` | control registry JSON | 4 CPU / 8G / 1h |
| `06_compare.sbatch` | patient vs control report | 4 CPU / 8G / 1h |
| `07_repertoire.sbatch` | genome-wide VNTR repertoire | 4 CPU / 16G / 2h |
| `08_detect_ceph_array.sbatch` | `detect` on the CEPH 1463 quartet | 16 CPU / 64G / 24h, x4 |
| `09_meiotic_instability.sbatch` | germline instability from CEPH duos | 4 CPU / 16G / 2h |
| `10_somatic_instability.sbatch` | per-read mosaicism across the cohort | 4 CPU / 32G / 4h |

## Germline instability: duo definition

`duos_ceph.tsv` defines the parent-offspring duos. For the CEPH 1463
quartet used here (NA12877 father, NA12878 mother, NA12879 daughter,
NA12882 son), this is a complete 2x2 design crossing parent sex with
offspring sex:

```
duo_id	parent_id	child_id	duo_type
CEPH_pat_79	NA12877	NA12879	father_daughter
CEPH_pat_82	NA12877	NA12882	father_son
CEPH_mat_79	NA12878	NA12879	mother_daughter
CEPH_mat_82	NA12878	NA12882	mother_son
```

Sex chromosomes are excluded by default (nearest-size allele matching
breaks down where either duo member is hemizygous).

Two outputs are produced. The **per-duo** summary, not the per-locus table,
is the correct unit of replication for comparing duo types: loci within one
duo share genetic background and sequencing run, so treating them as
independent observations would be pseudo-replication.

## Recommended order (first run, step by step)

1. `00_unit_tests.sbatch` — validates the package installation
2. `01_smoke_test_single_tool.sbatch` — validates `config.yaml`/paths (fast)
3. `02_detect_single_patient.sbatch` — validates the full pipeline on 1 patient
4. `03` / `04` / `08` — detection arrays (adjust `--array=1-N%K`)
5. `05_build_controls.sbatch` and `07_repertoire.sbatch` — after `03`
6. `06_compare.sbatch` — after `04` and `05`
7. `09_meiotic_instability.sbatch` — after `08`
8. `10_somatic_instability.sbatch` — after the relevant detect jobs (see below)

## Preparing the LongTR regions catalog

```bash
./scripts/prepare_longtr_bed.sh /mnt/references
```

Downloads the Adotto catalog (hg38) and converts it to LongTR's expected
format (`chrom, start 1-based, end, motif[,motif2], name`). Requires `gawk`.

## Chaining everything automatically

```bash
./scripts/submit_full_pipeline.sh
```

Submits 00 → 01 → (03, 04, 08 arrays in parallel) → 05, 07 → 06, 09, with
`--dependency=afterok:...` computed automatically and array sizes counted
dynamically from the TSV files.

`10_somatic_instability.sbatch` is deliberately **not** chained: it expects
all `{sample_id}/` subdirectories under a single `--detect-dir`, which
depends on how you lay out controls/patients/CEPH outputs. Either run it
once per cohort with a matching samples list, or symlink all sample
directories into one directory first, then submit it manually.
