"""
`meiotic-instability` and `somatic-instability` subcommands (paper Methods
2.10-2.11).

GERMLINE (MEIOTIC) instability: for each parent-offspring duo, at each
locus shared between the parent's and the child's merged VCF, each child
allele is matched to its nearest-sized parental allele (a standard
simplification in the STR mutation-rate literature -- it does NOT use true
parent-of-origin phasing). The size difference (child - matched parent
allele) is the instability estimate. Comparable only within a tool
(VAMOS/tandem-genotypes/LongTR size units are not interchangeable, see
merge.py), so matching and diffs are computed per tool.

SOMATIC (MITOTIC) instability: assessed from individual-read-level size
measurements at a single sample, independently of family structure, using
LongTR's FORMAT/ALLREADS field and tandem-genotypes' raw per-read length
list (both preserve per-read data; VAMOS does not, since it reports
per-haplotype consensus assemblies, and is therefore not used here). A read
is flagged as a candidate mosaic observation if it differs from every
called allele by at least one full repeat-motif unit -- a threshold chosen
to be robust to the ONT indel error rate in repetitive sequence. A locus is
only called mosaic if a minimum number of reads support the same off-allele
size, to avoid single-read sequencing errors being counted as instability.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path

import pandas as pd
import pysam

from str_toolkit import merge
from str_toolkit.annotate import DEFAULT_PROMOTER_WINDOW_BP, classify_location, classify_motif, load_exons, load_genes
from str_toolkit.utils import read_tsv_dicts

logger = logging.getLogger(__name__)

VALID_DUO_TYPES = {"mother_son", "mother_daughter", "father_son", "father_daughter"}

# LongTR's ALLREADS field may in principle include a sentinel bucket for
# reads that could not be confidently placed (observed as "-999" in
# HipSTR-family tutorials, HipSTR being LongTR's short-read ancestor).
# CHECKED against real LongTR output (multiple loci, including a large
# pathogenic expansion at C9orf72 with a 6721bp allele) -- no such sentinel
# value was observed. Kept as a defensive safeguard against future/other
# LongTR versions rather than removed: real bp-diff values seen so far
# range from -11 to 6721, so a -900 cutoff does not affect any genuine data.
ALLREADS_SENTINEL_THRESHOLD = -900


# ---------------------------------------------------------------------
# Germline (meiotic) instability
# ---------------------------------------------------------------------

def read_duos(path: str) -> list[dict]:
    """TSV with columns: duo_id, parent_id, child_id, duo_type."""
    duos = read_tsv_dicts(path, {"duo_id", "parent_id", "child_id", "duo_type"})
    invalid = sorted({d["duo_type"] for d in duos} - VALID_DUO_TYPES)
    if invalid:
        raise SystemExit(
            f"{path}: invalid duo_type value(s) {invalid}. "
            f"Expected one of {sorted(VALID_DUO_TYPES)}."
        )
    return duos


def match_transmitted_alleles(parent_sizes: list[float], child_sizes: list[float]) -> list[tuple[float, float, float]]:
    """
    For each child allele, finds its nearest-sized parental allele (putative
    transmitted allele). Returns a list of (child_allele, matched_parent_allele, diff).
    """
    if not parent_sizes or not child_sizes:
        return []
    matches = []
    for c in child_sizes:
        p = min(parent_sizes, key=lambda x: abs(x - c))
        matches.append((c, p, c - p))
    return matches


def compute_meiotic_instability(
    data_dir: Path,
    duos: list[dict],
    genes_bed: str,
    exons_bed: str,
    promoter_bp: int = DEFAULT_PROMOTER_WINDOW_BP,
) -> pd.DataFrame:
    data_dir = Path(data_dir)
    dict_genes = load_genes(genes_bed)
    dict_exons = load_exons(exons_bed)

    rows = []
    for duo in duos:
        duo_id, parent_id, child_id, duo_type = duo["duo_id"], duo["parent_id"], duo["child_id"], duo["duo_type"]

        parent_vcf = data_dir / parent_id / f"{parent_id}.merged.vcf"
        child_vcf = data_dir / child_id / f"{child_id}.merged.vcf"
        parent_loci = {f"{r['chrom']}_{r['pos']}_{r['motif']}": r for r in merge.parse_merged_vcf(parent_vcf)}
        child_loci = {f"{r['chrom']}_{r['pos']}_{r['motif']}": r for r in merge.parse_merged_vcf(child_vcf)}

        shared_loci = set(parent_loci) & set(child_loci)
        if not shared_loci:
            logger.warning("Duo %s: no shared loci between %s and %s", duo_id, parent_id, child_id)

        for locus_id in shared_loci:
            p_rec, c_rec = parent_loci[locus_id], child_loci[locus_id]

            p_by_tool: dict[str, list[float]] = defaultdict(list)
            for src, size in p_rec["sizes_by_source"].items():
                p_by_tool[merge.tool_family(src)].append(size)
            c_by_tool: dict[str, list[float]] = defaultdict(list)
            for src, size in c_rec["sizes_by_source"].items():
                c_by_tool[merge.tool_family(src)].append(size)

            loc_cat = classify_location(p_rec["chrom"], p_rec["pos"], dict_genes, dict_exons, promoter_bp)
            motif_cat = classify_motif(p_rec["motif"])

            for tool in sorted(set(p_by_tool) & set(c_by_tool)):
                for child_allele, parent_allele, diff in match_transmitted_alleles(p_by_tool[tool], c_by_tool[tool]):
                    rows.append({
                        "duo_id": duo_id,
                        "duo_type": duo_type,
                        "parent_id": parent_id,
                        "child_id": child_id,
                        "tool": tool,
                        "chrom": p_rec["chrom"],
                        "pos": p_rec["pos"],
                        "motif": p_rec["motif"],
                        "location_category": loc_cat,
                        "motif_category": motif_cat,
                        "parent_allele": parent_allele,
                        "child_allele": child_allele,
                        "diff": diff,
                    })

    return pd.DataFrame(rows)


def summarize_meiotic_instability_by_duo(df: pd.DataFrame) -> pd.DataFrame:
    """
    One row per (duo_id, tool, location_category, motif_category): median
    diff and number of loci. This -- not the raw per-locus table -- is the
    correct unit of replication for comparing duo types (mother-son vs.
    father-daughter, etc.): the thousands of loci within a single duo are
    not independent observations (shared genetic background, sequencing
    run), so treating them as such would understate the true uncertainty
    (pseudo-replication). Compare duo types using these per-duo summaries,
    e.g. with a Kruskal-Wallis test or a mixed model with duo as a random
    effect if working from the per-locus table directly.
    """
    if df.empty:
        return pd.DataFrame(columns=["duo_id", "duo_type", "tool", "location_category", "motif_category", "median_diff", "n_loci"])
    return (
        df.groupby(["duo_id", "duo_type", "tool", "location_category", "motif_category"])["diff"]
        .agg(median_diff="median", n_loci="count")
        .reset_index()
    )


# ---------------------------------------------------------------------
# Somatic (mitotic) instability
# ---------------------------------------------------------------------

def _parse_allreads(raw) -> list[tuple[float, int]]:
    """
    Parses LongTR's FORMAT/ALLREADS value, e.g. "-8|31;4|39" -> per-read
    bp-difference-from-reference buckets with their read counts. Excludes
    the sentinel bucket (see ALLREADS_SENTINEL_THRESHOLD).
    """
    if raw is None:
        return []
    if isinstance(raw, (tuple, list)):
        raw = ",".join(str(x) for x in raw if x is not None)
    pairs = []
    for chunk in str(raw).split(";"):
        if "|" not in chunk:
            continue
        bp_str, count_str = chunk.split("|", 1)
        try:
            bp, count = float(bp_str), int(count_str)
        except ValueError:
            continue
        if bp > ALLREADS_SENTINEL_THRESHOLD:
            pairs.append((bp, count))
    return pairs


def parse_longtr_for_somatic(vcf_path: Path) -> list[dict]:
    """
    Reads a LongTR VCF for somatic mosaicism analysis: called alleles
    (FORMAT/GB) plus the full per-read distribution (FORMAT/ALLREADS).
    """
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        return []
    records = []
    with pysam.VariantFile(str(vcf_path)) as vf:
        for record in vf:
            motif_field = record.info.get("MOTIF")
            if motif_field is None:
                continue
            motif = str(merge._first(motif_field)).split(",")[0]

            for sample_name, sample_data in record.samples.items():
                called_alleles = merge._parse_pipe_values(sample_data.get("GB"))
                allreads_raw = sample_data.get("ALLREADS")
                if not called_alleles or allreads_raw is None:
                    continue
                allreads = _parse_allreads(allreads_raw)
                if not allreads:
                    continue
                records.append({
                    "chrom": record.chrom,
                    "pos": record.pos,
                    "motif": motif,
                    "called_alleles": called_alleles,
                    "allreads": allreads,
                    "sample": sample_name,
                })
    return records


def detect_mosaicism(
    motif: str,
    called_alleles: list[float],
    allreads: list[tuple[float, int]],
    min_off_allele_reads: int = 3,
) -> dict:
    """
    A read is a candidate mosaic observation if it differs from every
    called allele by at least one full repeat-motif unit. A locus is only
    flagged mosaic if at least `min_off_allele_reads` reads support this.
    """
    motif_len = max(len(motif.strip()), 1)
    total_reads = sum(count for _, count in allreads)
    off_allele_reads = 0
    for bp_diff, count in allreads:
        if not called_alleles:
            continue
        nearest_dist = min(abs(bp_diff - allele) for allele in called_alleles)
        if nearest_dist >= motif_len:
            off_allele_reads += count

    return {
        "total_reads": total_reads,
        "off_allele_reads": off_allele_reads,
        "mosaic_fraction": (off_allele_reads / total_reads) if total_reads else 0.0,
        "is_mosaic": off_allele_reads >= min_off_allele_reads,
    }


def compute_somatic_instability(
    detect_dir: Path,
    sample_ids: list[str],
    genes_bed: str,
    exons_bed: str,
    promoter_bp: int = DEFAULT_PROMOTER_WINDOW_BP,
    min_off_allele_reads: int = 3,
) -> pd.DataFrame:
    detect_dir = Path(detect_dir)
    dict_genes = load_genes(genes_bed)
    dict_exons = load_exons(exons_bed)

    rows = []
    for sid in sample_ids:
        sample_dir = detect_dir / sid

        for rec in parse_longtr_for_somatic(sample_dir / f"{sid}.longtr.vcf.gz"):
            metrics = detect_mosaicism(rec["motif"], rec["called_alleles"], rec["allreads"], min_off_allele_reads)
            rows.append({
                "sample_id": sid, "tool": "longtr",
                "chrom": rec["chrom"], "pos": rec["pos"], "motif": rec["motif"],
                "location_category": classify_location(rec["chrom"], rec["pos"], dict_genes, dict_exons, promoter_bp),
                "motif_category": classify_motif(rec["motif"]),
                **metrics,
            })

        for row in merge._read_tandem_genotypes_rows(sample_dir / f"{sid}.tandem_genotypes.tsv"):
            allele1, allele2 = merge._split_two_alleles(row["values"])
            allreads = [(v, 1) for v in row["values"]]
            metrics = detect_mosaicism(row["motif"], [allele1, allele2], allreads, min_off_allele_reads)
            rows.append({
                "sample_id": sid, "tool": "tandem-genotypes",
                "chrom": row["chrom"], "pos": row["start"], "motif": row["motif"],
                "location_category": classify_location(row["chrom"], row["start"], dict_genes, dict_exons, promoter_bp),
                "motif_category": classify_motif(row["motif"]),
                **metrics,
            })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------
# CLI entry points
# ---------------------------------------------------------------------

def run_meiotic(args) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    duos = read_duos(args.duos)
    df = compute_meiotic_instability(args.data_dir, duos, args.genes_bed, args.exons_bed, args.promoter_bp)

    sep = "," if args.format == "csv" else "\t"
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep=sep, index=False)
    logger.info("Meiotic instability (per-locus) written: %s (%d rows)", output_path, len(df))

    if args.summary:
        summary_df = summarize_meiotic_instability_by_duo(df)
        summary_path = Path(args.summary)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(summary_path, sep=sep, index=False)
        logger.info("Per-duo summary written: %s (%d rows)", summary_path, len(summary_df))

    return 0


def run_somatic(args) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    sample_ids = [row["sample_id"] for row in read_tsv_dicts(args.samples_list, {"sample_id"})]
    df = compute_somatic_instability(
        args.detect_dir, sample_ids, args.genes_bed, args.exons_bed, args.promoter_bp, args.min_off_allele_reads
    )

    sep = "," if args.format == "csv" else "\t"
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep=sep, index=False)
    logger.info("Somatic instability written: %s (%d rows)", output_path, len(df))

    return 0
