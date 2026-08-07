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

- `--bam`: already-aligned BAM, used by VAMOS (clair3 + whatshap + vamos --contig).
- `--fastq`: raw merged fastq(.gz), used by TRGT/LongTR (minimap2 alignment,
  `map-ont` preset) and tandem-genotypes (last-train/lastal). TRGT and
  LongTR automatically reuse the same `.sorted.bam` if both run in the
  same job (no duplicate alignment).

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
registry and diffs, gene/feature annotation), and `repertoire` (genome-wide
VNTR repertoire, classified by genomic location and motif).

Not yet implemented: meiotic (parent-offspring duo) and somatic (per-read
mosaicism) instability analyses -- planned as a future subcommand.

Output formats **confirmed** on real files:
- TRGT: `INFO/MOTIFS`, `FORMAT/AL` (allele lengths in bp).
- tandem-genotypes: 8-column TSV (`chrom, start, end, motif, (ignored), .,
  per-read lengths, .`). The per-read length column is split into 2 groups
  (short/long allele) at the largest gap, and the median of each group is
  taken as the allele size (bp). Since the TRF bed often lists several
  overlapping candidate motifs for the same locus, these candidates are
  deduplicated, keeping the one covered by the most reads.
- LongTR: `INFO/MOTIF`, `FORMAT/GB` (bp difference from reference, per
  allele) — per the official gymrek-lab/LongTR README.

The micromamba environments referenced in `config.yaml` (clair3,
whatshap-env, vamos, trgt, longtr, last_env, tandem-env) must already exist
on the execution machine.

## License

MIT
