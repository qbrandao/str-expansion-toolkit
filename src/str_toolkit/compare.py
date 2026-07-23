"""
Sous-commande `compare`.

Compare les tailles STR des patients (VCF fusionné VAMOS+TRGT+tandem-genotypes,
produit par `detect`) au registre de contrôles (`build-controls`), un diff
PAR OUTIL disponible à chaque locus (les tailles ne sont pas comparables
entre outils : VAMOS = longueur en unités de motif, TRGT = longueur en bp,
tandem-genotypes = delta de copies vs référence).

Une ligne est gardée si au moins un outil dépasse le seuil. `n_tools_expanded`
compte combien d'outils confirment l'expansion (signal de confiance :
une expansion vue par 2-3 outils orthogonaux est plus fiable qu'une vue par
un seul). Tri décroissant sur `max_diff`.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pandas as pd

from str_toolkit import merge
from str_toolkit.annotate import annotate_locus, load_exons, load_genes

logger = logging.getLogger(__name__)

TOOLS = ("vamos", "trgt", "tandem-genotypes")

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
            logger.warning("VCF fusionné introuvable pour %s: %s", pid, merged_vcf)
            continue

        for record in merge.parse_merged_vcf(merged_vcf):
            motif = record["motif"]
            if triplet_only and len(motif) < 3:
                continue

            chrom, pos = record["chrom"], record["pos"]
            locus_id = f"{chrom}_{pos}_{motif}"
            control_entry = controls_registry.get(locus_id)
            if not control_entry:
                continue  # STR jamais observé chez les contrôles : pas de base de comparaison

            # Taille patient par outil = max sur les allèles/haplotypes disponibles pour cet outil
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
        raise SystemExit(f"Dossier introuvable : {patients_dir}")

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

    logger.info("Rapport écrit : %s (%d lignes)", args.output, len(df))
    return 0
