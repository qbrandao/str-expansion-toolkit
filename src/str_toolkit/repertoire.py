"""
`repertoire` subcommand.

Builds the genome-wide VNTR/STR repertoire (paper Methods 2.8-2.9, Results
3.2): every locus observed in a control cohort's merged VCFs, classified
into a mutually exclusive genomic location category and a motif-length
category, with per-tool control statistics (max size, number of
observations) attached.

Reuses `controls.collect_control_calls` for the per-locus, per-tool
registry, so `repertoire` and `build-controls` see the exact same set of
loci -- the repertoire is simply that registry enriched with the
location/motif classification and reshaped into a flat table for analysis
(counts by category, distribution plots, etc., done downstream in R/Python
notebooks -- this subcommand produces the input table, not the figures).
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from str_toolkit import controls
from str_toolkit.annotate import DEFAULT_PROMOTER_WINDOW_BP, classify_location, classify_motif, load_exons, load_genes

logger = logging.getLogger(__name__)


def build_repertoire(
    controls_dir: Path,
    genes_bed: str,
    exons_bed: str,
    sample_ids: list[str] | None = None,
    promoter_bp: int = DEFAULT_PROMOTER_WINDOW_BP,
) -> pd.DataFrame:
    registry = controls.collect_control_calls(controls_dir, sample_ids)
    dict_genes = load_genes(genes_bed)
    dict_exons = load_exons(exons_bed)

    rows = []
    for locus_id, entry in registry.items():
        chrom, pos, motif = entry["chrom"], entry["pos"], entry["motif"]
        row = {
            "locus_id": locus_id,
            "chrom": chrom,
            "pos": pos,
            "motif": motif,
            "location_category": classify_location(chrom, pos, dict_genes, dict_exons, promoter_bp),
            "motif_category": classify_motif(motif),
        }
        for tool, stats in entry["tools"].items():
            col = tool.replace("-", "_")
            row[f"{col}_max_size"] = stats["max_size"]
            row[f"{col}_n_observed"] = stats["n_observed"]
        rows.append(row)

    return pd.DataFrame(rows)


def summarize_repertoire(df: pd.DataFrame) -> pd.DataFrame:
    """Locus counts per (location_category, motif_category) cell -- paper Table 2."""
    return (
        df.groupby(["location_category", "motif_category"])
        .size()
        .reset_index(name="n_loci")
        .sort_values(["location_category", "motif_category"])
        .reset_index(drop=True)
    )


def run(args) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    controls_dir = Path(args.controls_dir)
    if not controls_dir.exists():
        raise SystemExit(f"Directory not found: {controls_dir}")

    sample_ids = None
    if args.samples_list:
        from str_toolkit.utils import read_samples_list

        sample_ids = [row["sample_id"] for row in read_samples_list(args.samples_list)]

    df = build_repertoire(
        controls_dir,
        genes_bed=args.genes_bed,
        exons_bed=args.exons_bed,
        sample_ids=sample_ids,
        promoter_bp=args.promoter_bp,
    )

    sep = "," if args.format == "csv" else "\t"
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep=sep, index=False)
    logger.info("Repertoire written: %s (%d loci)", output_path, len(df))

    if args.summary:
        summary_df = summarize_repertoire(df)
        summary_path = Path(args.summary)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(summary_path, sep=sep, index=False)
        logger.info("Summary written: %s (%d category cells)", summary_path, len(summary_df))

    return 0
