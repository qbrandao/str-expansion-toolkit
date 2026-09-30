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
| `11_validate_truthset.sbatch` | comparison against the Platinum Pedigree TR truthset | 4 CPU / 16G / 2h |
| `12_test_new_bam_longtr.sbatch` | fast check of a newly received BAM, LongTR alone | 8 CPU / 32G / 6h |
| `13_test_new_bam_full.sbatch` | full three-tool run on that same BAM | 16 CPU / 64G / 36h |

`inspect_bam.sh` is not a SLURM script: see "Testing a newly received BAM"
below.

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

## Testing a newly received BAM

A BAM from a sequencing provider or a collaborator is often the only input
available, with no FASTQ and no reliable description of how it was produced.
Three steps, in order.

**1. Triage, on the login node, well under a minute.**

```bash
./scripts/inspect_bam.sh /path/to/sample.bam /path/to/reference.fa > bam_report.txt
```

The second argument is optional and is the reference from `config.yaml`. Given
it, the report also compares the BAM contigs against the reference index, which
is the quickest way to catch a build mismatch. The report reads the header, the
index summary and a bounded sample of reads, so it is cheap enough for a login
node. Read sections 1, 2, 4 and 6 first:

| Finding | What it means |
|---|---|
| section 1, not coordinate-sorted | sort it before anything else |
| section 2, another build | pass `--realign`, which needs a FASTQ |
| section 2, contig naming differs from the catalogs | fix the catalogs, never the BAM |
| section 4, error rate at or above 3 percent | R9.4.1, so use an R9 clair3 model |
| section 4, mean read length under 500 bp | not long-read data, stop here |
| section 6, targeted | a genome-wide repertoire run is not meaningful |

**2. Fast single-tool check.**

```bash
sbatch scripts/12_test_new_bam_longtr.sbatch
```

LongTR alone on the BAM as it stands: no realignment, no clair3, no LAST. It
reports how many catalog loci were actually genotyped, which is what catches a
contig naming mismatch or insufficient depth before hours are spent on the slow
stages.

**3. Full run.**

```bash
sbatch scripts/13_test_new_bam_full.sbatch
```

All three default tools plus the merge, and it prints how many loci ended up
supported by one, two or three tools.

### A BAM is a sufficient input on its own

`detect` accepts `--bam`, `--fastq`, or both.

- With `--bam` alone, VAMOS, LongTR and TRGT read that BAM directly, with no
  realignment. tandem-genotypes needs reads rather than an alignment, so a
  FASTQ is extracted from the BAM once (primary alignments only, `MM`/`ML`
  methylation tags kept) and reused.
- With `--fastq` alone, the alignment is produced once with
  `minimap2 -ax map-ont -Y` and shared across the tools that need it.
- With both, the FASTQ wins wherever reads are needed, which is the better
  option: reads extracted from an aligned BAM are missing anything that failed
  to align, and LAST would never see them.
- `--realign` ignores `--bam` for the alignment-based tools and realigns from
  the FASTQ. Use it when the input BAM comes from another reference build,
  another aligner, or an unsuitable preset.

An input BAM with no index is never indexed in place, since it often sits in a
read-only or shared directory. It is symlinked into the sample output directory
and the symlink is indexed instead.

Indexing a BAM and extracting a FASTQ from it belong to no single tool, so
samtools for those steps comes from the `samtools_env` key in `config.yaml`, or
from the current `PATH` when that key is absent.

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
