"""
`compare` subcommand.

Compares patient STR sizes (merged VCF from VAMOS+TRGT+tandem-genotypes+LongTR,
produced by `detect`) to the control registry (`build-controls`), computing
a diff PER AVAILABLE TOOL at each locus (sizes are not comparable across
tools: VAMOS = length in motif-repeat units, TRGT = length in bp,
tandem-genotypes = length in bp derived from read clustering, LongTR = bp
difference from reference).

A row is kept if at least one tool exceeds the threshold. `n_tools_expanded`
counts how many tools confirm the expansion (a confidence signal: an
expansion seen by 2-3 orthogonal tools is more reliable than one seen by a
single tool). Sorted descending on `max_diff`.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pandas as pd

from str_toolkit import merge
from str_toolkit.annotate import annotate_locus, load_exons, load_genes

logger = logging.getLogger(__name__)

TOOLS = ("vamos", "trgt", "tandem-genotypes", "longtr")

OUTPUT_COLUMNS = (
    ["patient_id", "chrom", "pos", "motif", "gene", "feature"]
    + [f"{t.replace('-', '_')}_size" for t in TOOLS]
    + [f"{t.replace('-', '_')}_control_max" for t in TOOLS]
    + [f"{t.replace('-', '_')}_diff" for t in TOOLS]
    + ["n_tools_expanded", "max_diff"]
)


def build_comparison_table(
    patients_dir: Path,
    patient_ids: list[str],
    controls_registry: dict,
    genes_bed: str,
    exons_bed: str,
    threshold: int = 0,
    triplet_only: bool = False,
) -> pd.DataFrame:
    patients_dir = Path(patients_dir)
    dict_genes = load_genes(genes_bed)
    dict_exons = load_exons(exons_bed)

    rows = []
    for pid in patient_ids:
        merged_vcf = patients_dir / pid / f"{pid}.merged.vcf"
        if not merged_vcf.exists():
            logger.warning("Merged VCF not found for %s: %s", pid, merged_vcf)
            continue

        for record in merge.parse_merged_vcf(merged_vcf):
            motif = record["motif"]
            if triplet_only and len(motif) < 3:
                continue

            chrom, pos = record["chrom"], record["pos"]
            locus_id = f"{chrom}_{pos}_{motif}"
            control_entry = controls_registry.get(locus_id)
            if not control_entry:
                continue  # STR never observed in controls: no basis for comparison

            # Patient size per tool = max across the alleles/haplotypes available for that tool
            sizes_by_tool: dict[str, float] = {}
            for source, size in record["sizes_by_source"].items():
                tool = merge.tool_family(source)
                sizes_by_tool[tool] = max(sizes_by_tool.get(tool, size), size)

            diffs: dict[str, float] = {}
            row = {"patient_id": pid, "chrom": chrom, "pos": pos, "motif": motif}
            for tool in TOOLS:
                col = tool.replace("-", "_")
                patient_size = sizes_by_tool.get(tool)
                control_max = control_entry["tools"].get(tool, {}).get("max_size")
                row[f"{col}_size"] = patient_size
                row[f"{col}_control_max"] = control_max
                if patient_size is not None and control_max is not None:
                    diff = patient_size - control_max
                    row[f"{col}_diff"] = diff
                    diffs[tool] = diff
                else:
                    row[f"{col}_diff"] = None

            if not diffs:
                continue
            max_diff = max(diffs.values())
            if max_diff <= threshold:
                continue

            row["n_tools_expanded"] = sum(1 for d in diffs.values() if d > threshold)
            row["max_diff"] = max_diff

            genes, features = annotate_locus(chrom, pos, dict_genes, dict_exons)
            row["gene"] = genes
            row["feature"] = features

            rows.append(row)

    df = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    return df.sort_values("max_diff", ascending=False).reset_index(drop=True)


def run(args) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    with open(args.controls_json) as fh:
        controls_registry = json.load(fh)

    patients_dir = Path(args.patients_dir)
    if not patients_dir.exists():
        raise SystemExit(f"Directory not found: {patients_dir}")

    patient_ids = args.patients or sorted(p.name for p in patients_dir.iterdir() if p.is_dir())

    df = build_comparison_table(
        patients_dir,
        patient_ids,
        controls_registry,
        genes_bed=args.genes_bed,
        exons_bed=args.exons_bed,
        threshold=args.threshold,
        triplet_only=args.triplet_only,
    )

    sep = "," if args.format == "csv" else "\t"
    df.to_csv(args.output, sep=sep, index=False)

    logger.info("Report written: %s (%d rows)", args.output, len(df))
    return 0
