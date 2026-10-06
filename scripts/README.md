# SLURM scripts

Run from the repo root (`sbatch scripts/00_unit_tests.sbatch`, etc.), or
chain everything at once with `./scripts/submit_full_pipeline.sh`.

## Prerequisites (place at the repo root before running)

- an environment runner: `micromamba`, `mamba` or `conda`. Whichever is present
  is auto-detected, and `env_runner` in `config.yaml` overrides that with a name
  or a full path (an anaconda installation typically needs
  `env_runner: /home/<user>/anaconda3/bin/conda`)

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

`inspect_bam.sh`, `find_tools.sh` and `run_local.sh` are not SLURM scripts: see
"Testing a newly received BAM", "Finding the tools on a new machine" and
"Running without SLURM" below.

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

## Two annotation windows that are scientific choices

`repertoire` and both instability analyses take `--promoter-bp` and
`--subtelomere-bp`. Neither is a constant of nature, both belong in a methods
section, and both must stay identical across cases and controls.

`--subtelomere-bp` is measured from the chromosome ends, not from the
assembly's telomere annotation. That annotation is a 10 kb gap of N at each
end: no read aligns there, so no repeat can ever be called inside it. Using it
as the subtelomeric definition made the category structurally unreachable,
which a run on a real genome showed plainly, with zero subtelomeric loci and
the first callable locus of every chromosome sitting a few hundred base pairs
past the 10 kb mark. The default is 500 kb; `0` reverts to the gap and yields
nothing. On a real two-million-locus call set the window decides a meaningful
count:

| `--subtelomere-bp` | subtelomeric loci |
|---|---|
| 0 (the telomere gap) | 0 |
| 100 kb | 1 612 |
| 500 kb | 15 120 |

`--promoter-bp` is the window upstream of the TSS counted as `5prime_region`,
2 kb by default, where 1 kb and 5 kb are both defensible.

One more thing the categories depend on, and which no flag controls: how
inclusive the genes BED is. A file covering lncRNAs and pseudogenes pushes
loci from `intergenic_other` into `intronic`, so the intronic share describes
the annotation as much as the genome. State which gene set was used.

## Producing something a collaborator can read

A finished `detect` leaves a merged VCF of about two million rows whose
per-tool sizes are in three different units. It is not a report, and `compare`
cannot turn it into one yet, because `compare` answers "larger than in the
controls" and needs the registry that `build-controls` produces from a control
cohort.

`extract_loci_report.py` produces the one thing that is defensible from a single
sample: what each tool measured at a named set of genes.

```bash
./scripts/extract_loci_report.py \
  --merged-vcf results/local/1312_DFT/1312_DFT.merged.vcf \
  --genes-bed genes.bed.gz \
  --genes-file scripts/disease_genes_str.txt \
  --flank 5000 \
  -o 1312_DFT_report.tsv
```

`scripts/disease_genes_str.txt` lists 47 genes carrying a disease-associated
repeat, as a starting point rather than a curated clinical panel. Only gene
names are in it on purpose: repeat coordinates and pathogenic thresholds differ
between reference builds and between sources, so the coordinates come from the
genes BED and the thresholds stay the clinician's call.

The report carries a comment block stating that no locus is called expanded,
and gives both the raw value per tool and a `*_motif_units` conversion so the
columns can be read side by side:

| Tool | Raw unit | Converts how |
|---|---|---|
| VAMOS | motif-repeat units, per haplotype | already in motif units |
| LongTR | bp **difference** from the reference allele | divided by motif length, stays a difference |
| tandem-genotypes | bp from read-level clustering | divided by motif length |
| TRGT | absolute allele length in bp | divided by motif length |

`--min-tools` restricts to loci several tools agree on, and
`largest_motif_units_any_tool` sorts the table. That column mixes an absolute
size with a difference, so it ranks candidates and measures nothing.

## Finding the tools on a new machine

The environment names in `config.example.yaml` (`clair3`, `whatshap-env`,
`vamos`, `longtr`, `last_env`, `tandem-env`) are the ones from the machine this
project was developed on. Elsewhere they are named differently, or a tool is
installed outside any environment, and the preflight then reports a binary as
missing when it is only somewhere else.

```bash
./scripts/find_tools.sh
```

It lists the environments the runner knows about, says which one provides each
binary the pipeline needs, and prints the `config.yaml` environment block that
matches. Paste that block in, then rerun `run_local.sh --preflight`.

Two details in how it chooses:

- Between several environments holding the same binary it prefers the one named
  after it, since that is normally the environment installed for the tool rather
  than one that merely has it as a dependency. `whatshap` present in `clair3`,
  `nanocaller` and `whatshap` gives `whatshap`.
- For samtools and minimap2 it prefers the `PATH` and leaves the key out, since
  those are usually system installs and pinning them to an unrelated
  environment hides which version is really running.

A note on the preflight's own check, which is stricter than it looks. A conda
environment keeps the system directories on `PATH`, so `command -v samtools`
from inside an environment also finds `/usr/bin/samtools` and would report the
environment as providing a tool it does not contain. The preflight therefore
compares the resolved path against `CONDA_PREFIX` and distinguishes three
outcomes: `[ok]` for a binary inside the environment, `[!]` for one that
resolves outside it and will run but with whatever version the machine has, and
`[!!]` for one that is not found at all.

A missing binary does not necessarily block the run, since each tool needs only
its own:

| Running | Binaries needed |
|---|---|
| LongTR alone | `LongTR`, `minimap2`, `samtools` |
| tandem-genotypes alone | `lastal`, `last-train`, `last-split`, `tandem-genotypes`, `samtools` |
| VAMOS alone | `run_clair3.sh`, `whatshap`, `vamos`, `samtools` |
| TRGT (opt-in) | `trgt`, `minimap2`, `samtools` |

## Running without SLURM

Every `.sbatch` script reads `SLURM_CPUS_PER_TASK` and `SLURM_SUBMIT_DIR`, so
none of them works outside a cluster. `run_local.sh` does the same work with
plain shell on a single machine.

```bash
./scripts/run_local.sh --sample 1312_DFT --bam /path/sample.bam --background
```

It sizes the thread count from the machine (all cores minus one, capped at 16),
checks the inputs, the catalogs named in `config.yaml` and each micromamba
environment before starting anything slow, then runs under `setsid nohup` so the
run survives the terminal closing. `nice` and `ionice` keep a shared machine
usable while it runs.

| Option | Effect |
|---|---|
| `--sample ID` | sample identifier, required |
| `--bam FILE` | aligned, coordinate-sorted BAM |
| `--fastq FILE` | reads, preferred over extracting them from the BAM |
| `--config FILE` | default `config.yaml` |
| `--tools "A B"` | subset of tools, quoted |
| `--regions FILE` | reduced LongTR catalog, written into a config copy |
| `--threads N` | override the automatic count |
| `--outdir DIR` | default `results/local` |
| `--realign` | realign from the FASTQ rather than reusing the BAM |
| `--background` | detach and return immediately |
| `--status` | report which outputs exist and whether a run is active |
| `--preflight` | run the checks only |
| `--dry-run` | print the command that would run |

Recommended order on a machine without SLURM, which matters more there than on a
cluster because there is no queue to absorb a mistake:

```bash
# 1. check the environment before spending any time on it
./scripts/run_local.sh --sample 1312_DFT --bam /path/sample.bam --preflight

# 2. LongTR on a reduced catalog, minutes rather than hours
awk '$1=="chr9"' /path/longtr.regions.bed > quick_chr9.bed
./scripts/run_local.sh --sample 1312_DFT --bam /path/sample.bam \
  --tools longtr --regions quick_chr9.bed

# 3. LongTR on the full catalog, in the background
./scripts/run_local.sh --sample 1312_DFT --bam /path/sample.bam \
  --tools longtr --background

# 4. the slow tools, one at a time so a failure is easy to attribute
./scripts/run_local.sh --sample 1312_DFT --bam /path/sample.bam \
  --tools vamos --background
./scripts/run_local.sh --sample 1312_DFT --bam /path/sample.bam \
  --tools tandem-genotypes --background

# 5. once every tool has run, the same command with no --tools merges them
./scripts/run_local.sh --sample 1312_DFT --bam /path/sample.bam --background
```

### Memory, and why a killed stage used to go unnoticed

`lastal` memory grows with the thread count and with how much query sequence it
holds at once. At the defaults, on a whole-genome ONT sample, it has exhausted
RAM and swap and taken a machine down. Two config keys bound it:
`tandem_genotypes.last_threads` and `tandem_genotypes.lastal_batch_size`, the
latter mapping to `lastal -i`. On a 64 GB machine, `last_threads: 8` with
`lastal_batch_size: 1G` is a reasonable starting point.

Until this was fixed, such a kill was worse than a crash. The tool pipelines run
through `bash -c`, and without `pipefail` bash reports only the exit status of
the last command: a killed `lastal` in `lastal ... | last-split -m1 > out.maf`
left `last-split` exiting 0 on its short input, so the step looked successful
with a truncated MAF, tandem-genotypes built a TSV from partial alignments, and
the idempotence check then treated that TSV as finished for good. Every pipeline
now runs under `set -o pipefail`, and every produced file is written under a
temporary name and renamed only on a clean exit, so a truncated output can never
be mistaken for a finished step.

A run made before this fix is not automatically suspect. Check the log for
`exit code 0` and the MAF for a complete final alignment block, and compare the
last query name in the MAF against the last reads of the FASTQ: if lastal
consumed the whole input, the results stand. Only a run whose log shows a
non-zero exit, or a MAF ending mid-line, needs the MAF and the TSV deleted and
that step rerun.

### Reclaiming space after a finished run

A whole-genome run leaves roughly 330 GB behind for about 35 GB of actual
results. Once the outputs that matter are present, the intermediates can go:

| Keep | Why |
|---|---|
| `*_assembly.hap1.vcf`, `*_assembly.hap2.vcf` | VAMOS calls |
| `*.longtr.vcf.gz` | LongTR calls |
| `*.tandem_genotypes.tsv` | tandem-genotypes calls |
| `*.merged.vcf` | the merged locus set |
| `phased_merge_output.vcf.gz*` | the phasing, if VAMOS may be rerun |

| Delete | Size on a 95 GB BAM |
|---|---|
| `*_alignments.maf` | ~116 GB, only an input to the TSV |
| `*.reads.fastq.gz` | ~55 GB, re-extractable from the BAM |
| `*_haplotagged.bam`, `*_h1.bam`, `*_h2.bam` and their indexes | ~160 GB, only inputs to the haplotype VCFs |
| `tmp/`, `log/` | clair3 working files |

Deleting those is safe: each tool short circuits on its own final output, so a
later rerun of the same sample skips straight past the missing intermediates
rather than rebuilding them. The per-read analyses (`somatic-instability`) read
the LongTR VCF and the tandem-genotypes TSV, not the BAMs, so they are
unaffected.

### Disk space, which the preflight refuses to guess about

LAST writes its MAF as plain text carrying both aligned sequences, so for a
whole-genome ONT sample it commonly lands between 1.5 and 3 times the size of
the input BAM, and the whole file has to exist before tandem-genotypes reads it.
The preflight sizes this from the BAM and fails the run rather than letting it
die partway through `lastal` after hours of work. A 95 GB BAM needs roughly
190 to 330 GB free, counting the extracted FASTQ.

Three ways out when the estimate does not fit: point `--outdir` at a larger
filesystem, run the other tools locally and tandem-genotypes on the cluster, or
subsample the reads for a pipeline test rather than for publishable numbers. The
estimate is skipped when tandem-genotypes is not in `--tools`, so a LongTR or
VAMOS run is never blocked by it.

Follow a run with `--status`, or `tail -f logs/local_<sample>.log`. A second
launch for a sample already running is refused rather than allowed to collide.

Only one sample at a time per output directory. Two concurrent runs of the same
sample would fight over the same files, and the pid file guard only catches
those started through this script.

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
