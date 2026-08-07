"""
`build-controls` subcommand.

Reads the merged VCF (produced by `detect`, one per control sample:
{controls_dir}/{sample_id}/{sample_id}.merged.vcf) and builds a JSON
registry: for each locus, the maximum size observed SEPARATELY PER TOOL
(VAMOS / TRGT / tandem-genotypes / LongTR), since their units are not
comparable.

Output JSON format:
{
  "chr1_12345_AAAG": {
    "chrom": "chr1", "pos": 12345, "motif": "AAAG",
    "tools": {
      "vamos":  {"max_size": 42, "n_observed": 87},
      "trgt":   {"max_size": 38, "n_observed": 90},
      "tandem-genotypes": {"max_size": 3, "n_observed": 85}
    }
  },
  ...
}
"""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path

from str_toolkit import merge

logger = logging.getLogger(__name__)


def discover_control_samples(controls_dir: Path) -> list[str]:
    return sorted(p.name for p in controls_dir.iterdir() if p.is_dir())


def collect_control_calls(controls_dir: Path, sample_ids: list[str] | None = None) -> dict[str, dict]:
    controls_dir = Path(controls_dir)
    sample_ids = sample_ids or discover_control_samples(controls_dir)

    registry: dict[str, dict] = {}
    for sid in sample_ids:
        merged_vcf = controls_dir / sid / f"{sid}.merged.vcf"
        if not merged_vcf.exists():
            logger.warning("Merged VCF not found for %s: %s", sid, merged_vcf)
            continue

        for record in merge.parse_merged_vcf(merged_vcf):
            locus_id = f"{record['chrom']}_{record['pos']}_{record['motif']}"
            entry = registry.setdefault(
                locus_id,
                {"chrom": record["chrom"], "pos": record["pos"], "motif": record["motif"], "tools": {}},
            )
            for source, size in record["sizes_by_source"].items():
                tool = merge.tool_family(source)
                tool_entry = entry["tools"].setdefault(tool, {"max_size": size, "n_observed": 0})
                tool_entry["max_size"] = max(tool_entry["max_size"], size)
                tool_entry["n_observed"] += 1

    return registry


def run(args) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    controls_dir = Path(args.controls_dir)
    if not controls_dir.exists():
        raise SystemExit(f"Directory not found: {controls_dir}")

    sample_ids = None
    if args.samples_list:
        with open(args.samples_list, newline="") as fh:
            sample_ids = [row["sample_id"] for row in csv.DictReader(fh, delimiter="\t")]

    registry = collect_control_calls(controls_dir, sample_ids)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        json.dump(registry, fh, indent=2, sort_keys=True)

    logger.info("Control registry written: %s (%d loci)", output_path, len(registry))
    return 0
