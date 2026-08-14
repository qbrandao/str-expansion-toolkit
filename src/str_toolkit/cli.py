"""
Entry point for the `str-toolkit` CLI.

Six subcommands:
  1) detect               : runs VAMOS + tandem-genotypes + LongTR (TRGT
                             opt-in) for one or more samples, merges the
                             outputs into a final VCF.
  2) build-controls        : reads the merged VCFs of a control cohort and
                             builds a JSON registry {locus -> max size
                             observed per tool}.
  3) compare               : compares patient expansions to the control
                             registry, produces a sorted report.
  4) repertoire            : builds the genome-wide VNTR repertoire from a
                             control cohort, classified by genomic location
                             and motif.
  5) meiotic-instability    : germline instability from parent-offspring
                             duos (nearest-size allele matching).
  6) somatic-instability    : somatic (mosaic) instability from per-read
                             heterogeneity within single samples.

Usage:
  str-toolkit detect --sample p01 --bam p01.sorted.bam --fastq p01.merged.fastq.gz \
      --config config.yaml -o out/
  str-toolkit detect --samples-list patients.tsv --config config.yaml -o out/
  str-toolkit build-controls --controls-dir out/controls/ -o controls.json
  str-toolkit compare --patients-dir out/patients/ --controls-json controls.json \
      --genes-bed genes.bed.gz --exons-bed exons.bed.gz -o report.tsv
  str-toolkit repertoire --controls-dir out/controls/ \
      --genes-bed genes.bed.gz --exons-bed exons.bed.gz \
      -o repertoire.tsv --summary repertoire_summary.tsv
  str-toolkit meiotic-instability --duos duos.tsv --data-dir out/ \
      --genes-bed genes.bed.gz --exons-bed exons.bed.gz \
      -o meiotic.tsv --summary meiotic_by_duo.tsv
  str-toolkit somatic-instability --samples-list all_samples.tsv --detect-dir out/ \
      --genes-bed genes.bed.gz --exons-bed exons.bed.gz -o somatic.tsv
"""

from __future__ import annotations

import argparse
import sys

from str_toolkit import detect, controls, compare, repertoire, instability
from str_toolkit.annotate import DEFAULT_PROMOTER_WINDOW_BP


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="str-toolkit",
        description="Tandem repeat expansion detection and comparison (VAMOS / tandem-genotypes / LongTR / TRGT).",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # ---------------------------------------------------------------
    # 1) detect
    # ---------------------------------------------------------------
    p_detect = subparsers.add_parser(
        "detect",
        help="Run the detection tools for one or more samples and merge into a VCF.",
    )
    sample_group = p_detect.add_mutually_exclusive_group(required=True)
    sample_group.add_argument(
        "--sample", help="Sample ID for a single sample (see --bam/--fastq)."
    )
    sample_group.add_argument(
        "--samples-list",
        help="TSV file with columns: sample_id, bam_path, fastq_path (one sample per line).",
    )
    p_detect.add_argument("--bam", help="Already-aligned BAM (used by VAMOS/clair3).")
    p_detect.add_argument("--fastq", help="Raw merged fastq(.gz) (used by TRGT/LongTR and tandem-genotypes).")
    p_detect.add_argument(
        "--config",
        required=True,
        help="Path to config.yaml (reference paths, catalogs, micromamba environments).",
    )
    p_detect.add_argument("-o", "--outdir", required=True, help="Output directory.")
    p_detect.add_argument(
        "--tools",
        nargs="+",
        default=["vamos", "tandem-genotypes", "longtr"],
        choices=["vamos", "trgt", "tandem-genotypes", "longtr"],
        help="Subset of tools to run (default: vamos, tandem-genotypes, longtr -- all "
        "ONT-native). TRGT is NOT run by default: it has no official support for ONT "
        "data (designed for PacBio HiFi); it only runs if explicitly requested here, "
        "e.g. --tools vamos trgt tandem-genotypes longtr.",
    )
    p_detect.add_argument("--threads", type=int, default=4)
    p_detect.set_defaults(func=detect.run)

    # ---------------------------------------------------------------
    # 2) build-controls
    # ---------------------------------------------------------------
    p_controls = subparsers.add_parser(
        "build-controls",
        help="Build the reference JSON (max size per tool) from a control cohort's merged VCFs.",
    )
    p_controls.add_argument(
        "--controls-dir",
        required=True,
        help="Output directory of `detect` for the control cohort "
        "({controls-dir}/{sample_id}/{sample_id}.merged.vcf).",
    )
    p_controls.add_argument(
        "--samples-list",
        help="Optional: TSV (sample_id column) to restrict to the listed samples "
        "(default: all subdirectories of --controls-dir).",
    )
    p_controls.add_argument("-o", "--output", required=True, help="Path to the output JSON.")
    p_controls.set_defaults(func=controls.run)

    # ---------------------------------------------------------------
    # 3) compare
    # ---------------------------------------------------------------
    p_compare = subparsers.add_parser(
        "compare",
        help="Compare patient expansions to the control registry.",
    )
    p_compare.add_argument(
        "--patients-dir",
        required=True,
        help="Output directory of `detect` for the patients "
        "({patients-dir}/{sample_id}/{sample_id}.merged.vcf).",
    )
    p_compare.add_argument(
        "--patients",
        nargs="+",
        help="Optional: list of patient IDs to compare "
        "(default: all subdirectories of --patients-dir).",
    )
    p_compare.add_argument(
        "--controls-json", required=True, help="JSON generated by `build-controls`."
    )
    p_compare.add_argument(
        "--genes-bed", required=True, help="Gzipped BED: chrom/start/end/gene."
    )
    p_compare.add_argument(
        "--exons-bed", required=True, help="Gzipped BED: chrom/start/end/GENE_exonN (e.g. MANE Select)."
    )
    p_compare.add_argument(
        "-x", "--threshold", type=int, default=0,
        help="Additional threshold: only keep rows where diff_vs_controls > threshold (default: 0).",
    )
    p_compare.add_argument(
        "-t", "--triplet-only", action="store_true",
        help="Only keep STRs with a motif of 3 bp or longer.",
    )
    p_compare.add_argument(
        "-o", "--output", required=True, help="Output file (TSV/CSV)."
    )
    p_compare.add_argument(
        "--format", choices=["tsv", "csv"], default="tsv", help="Output report format."
    )
    p_compare.set_defaults(func=compare.run)

    # ---------------------------------------------------------------
    # 4) repertoire
    # ---------------------------------------------------------------
    p_repertoire = subparsers.add_parser(
        "repertoire",
        help="Build the genome-wide VNTR repertoire (location/motif classification) from a control cohort.",
    )
    p_repertoire.add_argument(
        "--controls-dir",
        required=True,
        help="Output directory of `detect` for the control cohort "
        "({controls-dir}/{sample_id}/{sample_id}.merged.vcf).",
    )
    p_repertoire.add_argument(
        "--samples-list",
        help="Optional: TSV (sample_id column) to restrict to the listed samples "
        "(default: all subdirectories of --controls-dir).",
    )
    p_repertoire.add_argument(
        "--genes-bed", required=True, help="Gzipped BED: chrom/start/end/gene."
    )
    p_repertoire.add_argument(
        "--exons-bed", required=True, help="Gzipped BED: chrom/start/end/GENE_exonN (e.g. MANE Select)."
    )
    p_repertoire.add_argument(
        "--promoter-bp", type=int, default=DEFAULT_PROMOTER_WINDOW_BP,
        help=f"Promoter window (bp) upstream of the TSS counted as '5prime_region' (default: {DEFAULT_PROMOTER_WINDOW_BP}).",
    )
    p_repertoire.add_argument(
        "-o", "--output", required=True, help="Output file: one row per locus (TSV/CSV)."
    )
    p_repertoire.add_argument(
        "--summary",
        help="Optional: also write a location x motif category count summary to this path (paper Table 2).",
    )
    p_repertoire.add_argument(
        "--format", choices=["tsv", "csv"], default="tsv", help="Output format."
    )
    p_repertoire.set_defaults(func=repertoire.run)

    # ---------------------------------------------------------------
    # 5) meiotic-instability
    # ---------------------------------------------------------------
    p_meiotic = subparsers.add_parser(
        "meiotic-instability",
        help="Germline instability from parent-offspring duos (nearest-size allele matching).",
    )
    p_meiotic.add_argument(
        "--duos", required=True,
        help="TSV with columns: duo_id, parent_id, child_id, duo_type "
        "(duo_type in mother_son/mother_daughter/father_son/father_daughter).",
    )
    p_meiotic.add_argument(
        "--data-dir", required=True,
        help="Output directory of `detect` containing {sample_id}/{sample_id}.merged.vcf "
        "for every parent_id/child_id referenced in --duos.",
    )
    p_meiotic.add_argument("--genes-bed", required=True, help="Gzipped BED: chrom/start/end/gene.")
    p_meiotic.add_argument("--exons-bed", required=True, help="Gzipped BED: chrom/start/end/GENE_exonN.")
    p_meiotic.add_argument(
        "--promoter-bp", type=int, default=DEFAULT_PROMOTER_WINDOW_BP,
        help=f"Promoter window (bp) upstream of the TSS (default: {DEFAULT_PROMOTER_WINDOW_BP}).",
    )
    p_meiotic.add_argument("-o", "--output", required=True, help="Per-locus long-format output (TSV/CSV).")
    p_meiotic.add_argument(
        "--summary",
        help="Optional: also write a per-duo summary (median diff, n_loci) to this path -- "
        "the correct unit of replication for comparing duo types, see instability.py docstring.",
    )
    p_meiotic.add_argument("--format", choices=["tsv", "csv"], default="tsv", help="Output format.")
    sex_chrom_group = p_meiotic.add_mutually_exclusive_group()
    sex_chrom_group.add_argument(
        "--exclude-sex-chromosomes", dest="include_sex_chromosomes", action="store_false",
        help="Exclude X/Y loci (DEFAULT). Nearest-size allele matching assumes two comparable "
        "alleles per locus in both duo members, which fails wherever either is hemizygous "
        "(e.g. a son on X and Y), producing ploidy artefacts rather than instability.",
    )
    sex_chrom_group.add_argument(
        "--include-sex-chromosomes", dest="include_sex_chromosomes", action="store_true",
        help="Keep X/Y loci in the analysis. Only appropriate if hemizygosity is handled downstream.",
    )
    p_meiotic.set_defaults(include_sex_chromosomes=False)
    p_meiotic.set_defaults(func=instability.run_meiotic)

    # ---------------------------------------------------------------
    # 6) somatic-instability
    # ---------------------------------------------------------------
    p_somatic = subparsers.add_parser(
        "somatic-instability",
        help="Somatic (mosaic) instability from per-read heterogeneity within single samples.",
    )
    p_somatic.add_argument(
        "--samples-list", required=True,
        help="TSV with a sample_id column (any cohort: controls, patients, duo individuals -- "
        "no family structure required).",
    )
    p_somatic.add_argument(
        "--detect-dir", required=True,
        help="Output directory of `detect` containing {sample_id}/{sample_id}.longtr.vcf.gz "
        "and {sample_id}/{sample_id}.tandem_genotypes.tsv (the raw per-tool outputs, NOT the "
        "merged VCF, since read-level detail is lost during merging). VAMOS is not used here "
        "as it reports per-haplotype consensus rather than individual-read measurements.",
    )
    p_somatic.add_argument("--genes-bed", required=True, help="Gzipped BED: chrom/start/end/gene.")
    p_somatic.add_argument("--exons-bed", required=True, help="Gzipped BED: chrom/start/end/GENE_exonN.")
    p_somatic.add_argument(
        "--promoter-bp", type=int, default=DEFAULT_PROMOTER_WINDOW_BP,
        help=f"Promoter window (bp) upstream of the TSS (default: {DEFAULT_PROMOTER_WINDOW_BP}).",
    )
    p_somatic.add_argument(
        "--min-off-allele-reads", type=int, default=3,
        help="Minimum number of reads supporting the same off-allele size to call a locus "
        "mosaic (default: 3, to avoid single-read sequencing errors being counted as instability).",
    )
    p_somatic.add_argument("-o", "--output", required=True, help="Per-locus output (TSV/CSV).")
    p_somatic.add_argument("--format", choices=["tsv", "csv"], default="tsv", help="Output format.")
    p_somatic.set_defaults(func=instability.run_somatic)

    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args) or 0


if __name__ == "__main__":
    sys.exit(main())
