# str-expansion-toolkit

CLI for detecting and analyzing short tandem repeat (STR) expansions from
ONT long-read data, combining several orthogonal tools: **VAMOS**,
**tandem-genotypes**, **LongTR** (default), and **TRGT** (optional, see
note below).

## ⚠️ Important note on TRGT and ONT data

TRGT is designed for PacBio HiFi reads and **has no official support for
ONT data** (see Aliyev et al. 2026, bioRxiv, which explicitly excludes TRGT
from an ONT tandem repeat genotyping benchmark for this reason — it does
not accept ONT data as an intended use case). As a result, **TRGT is not
run by default** by `detect`. It remains available via `--tools ... trgt
...` for comparison/exploration, but any result produced this way must be
documented as an off-label use in any resulting publication.

The third default tool is **LongTR**, a long-read adaptation of HipSTR
designed for both PacBio HiFi **and** ONT reads.

## Installation

```bash
git clone https://github.com/<you>/str-expansion-toolkit.git
cd str-expansion-toolkit
pip install -e .
```

## Configuration

All installation-specific paths (reference genome, VAMOS motif catalog,
minimap2 index for TRGT/LongTR, LAST index for tandem-genotypes, clair3
model, micromamba environment names) are centralized in `config.yaml`,
based on the `config.example.yaml` template:

```bash
cp config.example.yaml config.yaml
# edit config.yaml with your own paths
```

Two optional keys cover how the machine is set up:

- `samtools_env` names the environment providing samtools for the two steps that
  belong to no single tool: indexing an input BAM that has no index, and
  extracting a FASTQ from a BAM when none is supplied. Left out, samtools is
  taken from the current `PATH`.
- `env_runner` is how the tool environments are entered: `micromamba`, `mamba`,
  `conda`, or a full path to one of them such as
  `/home/user/anaconda3/bin/conda`. Left out, the first of the three found on
  `PATH` is used. conda is called with `--no-capture-output`, since it otherwise
  buffers a long tool's progress away.
- `align_env` provides minimap2 and samtools for the shared alignment step,
  which only runs when no BAM is given or with `--realign`. Left out, the
  environment of the tool that asked for the alignment is used, which then has
  to contain minimap2 itself.

Each tool section also accepts the binary itself, as `bin`, `bin_vamos`,
`bin_clair3`, `bin_whatshap`, `bin_lastal` and so on. A bare name is resolved
inside the environment, which is the default. An absolute path is invoked as is,
which is what a tool compiled by hand outside any environment needs, vamos and
LongTR being commonly built that way:

```yaml
vamos:
  env_vamos: ""                                # no environment
  bin_vamos: /home/user/vamos/src/vamos        # hand-built
```

Setting the env to `""` runs the binary on the current `PATH`. Keeping the env
name alongside an absolute path also works, and is worth doing when the binary
needs that environment's shared libraries.

`scripts/find_tools.sh` works all of this out for a given machine, including
searching for hand-built binaries, and prints the block to paste in.

## Usage

### 1) Detect expansions in one or more patients

```bash
# A single patient (default tools: vamos, tandem-genotypes, longtr)
str-toolkit detect \
  --sample patient01 \
  --bam patient01.sorted.bam \
  --fastq patient01_Guppy_4.0.11_prom.merged.fastq.gz \
  --config config.yaml \
  -o results/patients/

# A list of patients (TSV: sample_id, bam_path, fastq_path)
str-toolkit detect \
  --samples-list patients.tsv \
  --config config.yaml \
  -o results/patients/

# Explicitly include TRGT (off-label on ONT data, only do this knowingly
# -- see the note at the top of this file)
str-toolkit detect \
  --sample patient01 --bam ... --fastq ... --config config.yaml \
  --tools vamos tandem-genotypes longtr trgt \
  -o results/patients/
```

#### Inputs: a BAM, a FASTQ, or both

Either input is sufficient on its own.

- `--bam`: an already-aligned, coordinate-sorted BAM. Read directly by VAMOS
  (clair3 + whatshap + `vamos --contig`), and reused as is by LongTR and TRGT,
  with no realignment. This assumes it was aligned to the `reference` given in
  `config.yaml`.
- `--fastq`: raw merged fastq(.gz). Needed by tandem-genotypes
  (last-train/lastal), which realigns the reads itself, and used to produce the
  alignment when no BAM is given. LongTR and TRGT share that single
  `.sorted.bam` when both run in the same job, so nothing is aligned twice.
- **BAM alone**: the reads tandem-genotypes needs are extracted from the BAM
  once (primary alignments only, `MM`/`ML` methylation tags kept) and reused.
- **Both**: the FASTQ wins wherever reads are needed, which is the better
  option. Reads extracted from an aligned BAM are missing anything that failed
  to align, and LAST would never see them.
- `--realign`: ignores `--bam` for the alignment-based tools and realigns from
  the FASTQ with `minimap2 -ax map-ont -Y`. Use it when the input BAM comes from
  another reference build, another aligner, or an unsuitable preset.

An input BAM with no index is never indexed in place, since it often sits in a
read-only or shared directory: it is symlinked into the sample output directory
and the symlink is indexed instead.

Before running any of this on a BAM whose provenance is not fully known, triage
it first:

```bash
./scripts/inspect_bam.sh /path/to/sample.bam /path/to/reference.fa > bam_report.txt
```

The report settles whether the file is long-read data at all, which reference
build and contig naming it uses, whether the chemistry is R9.4.1 or R10.4.1
(from the header, and independently from the observed error rate), the depth,
and whether it is whole genome or targeted. It reads only the header, the index
summary and a bounded sample of reads, so it is cheap enough to run on a login
node. See `scripts/README.md` for how to act on each finding, and
`scripts/12_test_new_bam_longtr.sbatch` then
`scripts/13_test_new_bam_full.sbatch` for the two test runs that follow.

Runs the selected tools for each sample (with fault tolerance: steps whose
output already exists are skipped), then merges the outputs into a summary
VCF: `results/patients/<sample_id>/<sample_id>.merged.vcf`.

**Merging (`str_toolkit/merge.py`):** the tools do not anchor their
coordinates the same way for a given biological locus. Merging therefore
groups calls by **tolerant interval** (`window`, 25 bp by default) **and
canonical motif** (circular rotation: `AAAG` ≡ `AAGA` ≡ `GAAA`), rather than
by strict positional equality. Sizes remain **source-specific** in the VCF
(`SIZES=vamos_hap1:42,longtr_allele1:12,...`) since their units are not
comparable:
- VAMOS: length in motif-repeat units
- TRGT: length in bp (absolute allele length)
- LongTR: bp difference from the reference (delta, not an absolute length)
- tandem-genotypes: length in bp derived from read-level length clustering

### 2) Build the control registry

```bash
str-toolkit build-controls \
  --controls-dir results/controls/ \
  -o results/controls.json
```

Reads each control sample's merged VCF
(`{controls-dir}/{sample_id}/{sample_id}.merged.vcf`, produced by `detect`)
and builds, for each locus, the maximum observed size **separately per
tool**:

```json
{
  "chr1_12345_AAAG": {
    "chrom": "chr1", "pos": 12345, "motif": "AAAG",
    "tools": {
      "vamos": {"max_size": 42, "n_observed": 87},
      "tandem-genotypes": {"max_size": 3, "n_observed": 85},
      "longtr": {"max_size": 10, "n_observed": 88}
    }
  }
}
```

Optional: `--samples-list` to restrict to the listed samples (otherwise all
subdirectories of `--controls-dir` are used).

### 3) Compare patients to controls

```bash
str-toolkit compare \
  --patients-dir results/patients/ \
  --controls-json results/controls.json \
  --genes-bed genes.bed.gz \
  --exons-bed MANE_Select_exons.bed.gz \
  -x 0 \
  -o report.tsv
```

For each patient locus present in the control registry, computes a diff
**per available tool** (`vamos_diff`, `trgt_diff`, `tandem_genotypes_diff`,
`longtr_diff`), keeps the row if at least one tool exceeds the threshold
(`-x/--threshold`), and sorts descending on `max_diff` (the largest diff
among the tools available for that locus).

Output columns:
`patient_id, chrom, pos, motif, gene, feature, vamos_size, vamos_control_max,
vamos_diff, trgt_size, trgt_control_max, trgt_diff, tandem_genotypes_size,
tandem_genotypes_control_max, tandem_genotypes_diff, longtr_size,
longtr_control_max, longtr_diff, n_tools_expanded, max_diff`

`n_tools_expanded` counts how many orthogonal tools confirm the expansion at
that locus — a useful confidence signal (an expansion seen by several
independent tools is more reliable than one seen by a single tool).

Options: `-t/--triplet-only` (motifs ≥ 3 bp only), `--patients` (restrict to
a list of IDs), `--format csv`.

### 4) Build the genome-wide VNTR repertoire

```bash
str-toolkit repertoire \
  --controls-dir results/controls/ \
  --genes-bed genes.bed.gz \
  --exons-bed MANE_Select_exons.bed.gz \
  --promoter-bp 2000 \
  -o repertoire.tsv \
  --summary repertoire_summary.tsv
```

Reads the same control cohort registry as `build-controls`, and classifies
every locus into ONE mutually exclusive genomic location category
(`subtelomeric`, `paracentromeric`, `5prime_region` [promoter window + 5'
UTR], `exonic`, `intronic`, `intergenic_other`) and one motif-length
category (`mononucleotide` ... `hexanucleotide_or_longer`), using
`str_toolkit/annotate.py`.

- `repertoire.tsv`: one row per locus, with `location_category`,
  `motif_category`, and per-tool `{tool}_max_size`/`{tool}_n_observed`
  columns.
- `--summary`: locus counts per (location, motif) category cell (paper
  Table 2).

`str_toolkit/annotate.py` also provides `annotate_locus(chrom, pos,
dict_genes, dict_exons)`, used internally by `compare` for per-patient
gene/feature reporting (exon, UTR, intronic, intergenic, centromere,
telomere) -- distinct from the mutually-exclusive `classify_location` used
by `repertoire`.

### 5) Germline (meiotic) instability from parent-offspring duos

```bash
str-toolkit meiotic-instability \
  --duos duos.tsv \
  --data-dir results/ \
  --genes-bed genes.bed.gz \
  --exons-bed MANE_Select_exons.bed.gz \
  -o meiotic_instability.tsv \
  --summary meiotic_instability_by_duo.tsv
```

`duos.tsv` columns: `duo_id, parent_id, child_id, duo_type` (`duo_type` in
`mother_son`/`mother_daughter`/`father_son`/`father_daughter`). Both
`parent_id` and `child_id` must have already been run through `detect`
(`{data-dir}/{id}/{id}.merged.vcf`).

**Sex chromosomes are excluded by default.** Nearest-size allele matching
assumes both duo members carry two comparable alleles at a locus, which
fails wherever either is hemizygous — a son is hemizygous for X and carries
a single Y — so X/Y loci would otherwise yield ploidy artefacts rather than
instability estimates. The number of excluded loci is logged per duo. Pass
`--include-sex-chromosomes` to keep them, only if hemizygosity is handled
downstream.

For each duo and each locus shared between parent and child, each child
allele is matched to its nearest-sized parental allele (a standard
simplification -- not true parent-of-origin phasing), and the size
difference is recorded, per available tool (VAMOS/tandem-genotypes/LongTR;
never merged across tools, see [Genomic and motif classification](#4-build-the-genome-wide-vntr-repertoire)
above on size units).

`-o` writes one row per (duo, locus, tool). `--summary` additionally writes
one row per (duo, tool, location, motif category) with the median diff and
locus count -- **this is the correct table to use when comparing duo
types** (mother-son vs. father-daughter, etc.): the thousands of loci
within a single duo share genetic background and sequencing run, so
treating each locus as an independent observation across duos would
understate uncertainty (pseudo-replication).

### 6) Somatic (mitotic) instability from per-read heterogeneity

```bash
str-toolkit somatic-instability \
  --samples-list all_samples.tsv \
  --detect-dir results/ \
  --genes-bed genes.bed.gz \
  --exons-bed MANE_Select_exons.bed.gz \
  --min-off-allele-reads 3 \
  -o somatic_instability.tsv
```

Unlike meiotic instability, this needs no family structure -- any sample
run through `detect` can contribute (controls, patients, duo individuals).
It reads the tools' **raw, per-sample outputs** directly
(`{detect-dir}/{id}/{id}.longtr.vcf.gz` and
`{detect-dir}/{id}/{id}.tandem_genotypes.tsv`), not the merged VCF, since
individual-read detail is collapsed away during merging. VAMOS is not used
here: it reports per-haplotype consensus assemblies, not per-read
measurements.

At each locus, a read is flagged as a candidate mosaic observation if it
differs from every called allele by at least one full repeat-motif unit
(robust to the ONT indel error rate in repetitive sequence). A locus is
only called mosaic if `--min-off-allele-reads` (default 3) reads support
the same off-allele size, to avoid single-read sequencing errors being
counted as instability.

The LongTR `FORMAT/ALLREADS` parser excludes a sentinel bucket inferred
from HipSTR-family tutorials (`bp_diff <= -900`, LongTR's short-read
ancestor) for reads that could not be confidently placed. **Checked
against real LongTR output (multiple loci) and no such sentinel value was
observed** -- kept as a defensive safeguard only; real bp-diff values seen
so far (-11 to 6721) are unaffected by this cutoff.

### 7) Validate against the Platinum Pedigree TR truthset

```bash
# Download once
aws s3 cp --no-sign-request \
  s3://platinum-pedigree-data/variants/tr_truthset/GRCh38/ceph_1463_tandem_repeats.oa.vcf.gz \
  references/

# Confirm which INFO field carries the motif before trusting any matching
str-toolkit validate-truthset --truthset references/ceph_1463_tandem_repeats.oa.vcf.gz \
  --inspect -o /dev/null

str-toolkit validate-truthset \
  --truthset references/ceph_1463_tandem_repeats.oa.vcf.gz \
  --instability results/meiotic_instability.tsv \
  -o results/meiotic_instability.truthset_annotated.tsv \
  --summary results/truthset_false_positive_rate.tsv
```

The Platinum Pedigree truthset (Kronenberg et al., Nat Methods 2025) is not
a list of de novo mutations. It is a set of tandem repeat genotypes retained
because they segregate consistently with Mendelian inheritance across CEPH
1463. At those loci transmission is therefore known to be consistent, so a
parent-to-child size difference called by `meiotic-instability` is most
likely a false positive. The per-tool `apparent_fp_rate` in the summary is
an **upper bound** on the false-positive rate of the nearest-size
transmitted-allele assignment: a locus may still be genuinely unstable and
pass the pedigree filters.

Loci are matched to the truthset with the same interval tolerance used at
the merge step (`--window`, 25 bp by default), since the truthset does not
share the coordinate convention of any single caller. Motifs are compared
after canonicalization; if the truthset carries no usable motif field the
matching falls back to position only and logs a warning.

## Running without SLURM

`scripts/run_local.sh` runs `detect` on a single machine with no scheduler,
under `nohup`, with a preflight of the environments and catalogs and a
`--status` view of progress. See `scripts/README.md`, section "Running without
SLURM".

```bash
./scripts/run_local.sh --sample 1312_DFT --bam /path/sample.bam --preflight
./scripts/run_local.sh --sample 1312_DFT --bam /path/sample.bam --background
./scripts/run_local.sh --sample 1312_DFT --status
```

## Sample file format

TSV file with a header, used by `--samples-list` in `detect` (`bam_path`
and/or `fastq_path` depending on which tools are used):

```
sample_id	bam_path	fastq_path
patient01	/data/bam/patient01.sorted.bam	/data/fastq/patient01.merged.fastq.gz
patient02	/data/bam/patient02.sorted.bam	/data/fastq/patient02.merged.fastq.gz
```

## Project status

End-to-end functional: `detect` (VAMOS/tandem-genotypes/LongTR by default,
TRGT opt-in, + multi-tool merging), `build-controls`, `compare` (per-tool
registry and diffs, gene/feature annotation), `repertoire` (genome-wide
VNTR repertoire, classified by genomic location and motif),
`meiotic-instability` (parent-offspring duos, nearest-size allele matching),
`somatic-instability` (per-read mosaicism detection via LongTR ALLREADS and
tandem-genotypes raw read lengths), and `validate-truthset` (empirical
false-positive rate against the Platinum Pedigree TR truthset).

All four tools' output formats are confirmed against real data (see below),
including LongTR's GB/ALLREADS pipe-separated encoding and its (absent, but
defensively handled) sentinel value.

Output formats **confirmed** on real files:
- VAMOS: `INFO/RU` (comma-separated candidate motifs, first one used),
  `INFO/LEN_H1` (size in motif-repeat units) -- **confirmed on a real
  VAMOS line for C9orf72 (chr9:27573415, LEN_H1=1123 units)**, cross-checked
  against the LongTR regression for the same gene (GB=6721bp / 6bp motif ≈
  1120 units -- consistent).
- TRGT: `INFO/MOTIFS`, `FORMAT/AL` (allele lengths in bp).
- tandem-genotypes: 8-column TSV (`chrom, start, end, motif, (ignored), .,
  forward-strand per-read values, reverse-strand per-read values`).
  Columns 7 and 8 (index 6/7) are **confirmed** to both carry per-read
  copy-number-change values when reverse-strand reads exist (real C9orf72
  line, chr9:27573484) -- both are pooled into a single per-read list, not
  just column 7 alone (an earlier bug silently dropped reverse-strand
  reads whenever present). The pooled list is split into 2 groups
  (short/long allele) at the largest gap, and the median of each group is
  taken as the allele size for this tool. Since the TRF bed often lists
  several overlapping candidate motifs for the same locus, these
  candidates are deduplicated, keeping the one covered by the most reads.
- LongTR: `INFO/MOTIF`, `FORMAT/GB` and `FORMAT/ALLREADS` -- **confirmed on
  two independent real LongTR lines (C9orf72, chr9:27573455; and a second
  locus, chr3:194943234)**. Both fields pack their per-allele/per-read
  values into a single `|`-joined string (e.g. `GB=-6|6721`), NOT the
  comma-separated array the VCF spec would suggest for a multi-valued
  FORMAT field. Naive iteration over the raw pysam value would silently
  iterate over string characters instead of alleles -- always parse GB via
  `merge._parse_pipe_values`, never index into it directly.

The micromamba environments referenced in `config.yaml` (clair3,
whatshap-env, vamos, trgt, longtr, last_env, tandem-env) must already exist
on the execution machine.

## License

MIT
