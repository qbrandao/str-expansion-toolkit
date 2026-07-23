"""
Point d'entrée du CLI `str-toolkit`.

Trois sous-commandes :
  1) detect         : lance VAMOS + TRGT + tandem-genotypes pour un ou
                       plusieurs patients, fusionne les sorties en un VCF final.
  2) build-controls  : lance la détection sur une cohorte de contrôles et
                       génère un JSON {STR -> taille max observée, ...}.
  3) compare         : compare les expansions des patients au JSON de
                       contrôles, produit un tableau trié.

Usage :
  str-toolkit detect --sample p01 --bam p01.sorted.bam --fastq p01.merged.fastq.gz \
      --config config.yaml -o out/
  str-toolkit detect --samples-list patients.tsv --config config.yaml -o out/
  str-toolkit build-controls --samples-list controls.tsv --config config.yaml -o controls.json
  str-toolkit compare --patients-vcf out/*/*.merged.vcf --controls-json controls.json -o report.tsv
"""

from __future__ import annotations

import argparse
import sys

from str_toolkit import detect, controls, compare


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="str-toolkit",
        description="Détection et comparaison d'expansions STR (VAMOS / TRGT / tandem-genotypes).",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # ---------------------------------------------------------------
    # 1) detect
    # ---------------------------------------------------------------
    p_detect = subparsers.add_parser(
        "detect",
        help="Lance les 3 outils sur un ou plusieurs patients et fusionne en un VCF.",
    )
    sample_group = p_detect.add_mutually_exclusive_group(required=True)
    sample_group.add_argument(
        "--sample", help="Identifiant + chemin BAM/CRAM d'un seul patient (voir --bam)."
    )
    sample_group.add_argument(
        "--samples-list",
        help="Fichier TSV avec colonnes: sample_id, bam_path (un patient par ligne).",
    )
    p_detect.add_argument("--bam", help="BAM déjà aligné (utilisé par VAMOS/clair3).")
    p_detect.add_argument("--fastq", help="Fastq(.gz) brut fusionné (utilisé par TRGT et tandem-genotypes).")
    p_detect.add_argument(
        "--config",
        required=True,
        help="Fichier config.yaml (chemins refs, catalogues, environnements micromamba).",
    )
    p_detect.add_argument("-o", "--outdir", required=True, help="Dossier de sortie.")
    p_detect.add_argument(
        "--tools",
        nargs="+",
        default=["vamos", "trgt", "tandem-genotypes"],
        choices=["vamos", "trgt", "tandem-genotypes"],
        help="Sous-ensemble d'outils à lancer (par défaut: les 3).",
    )
    p_detect.add_argument("--threads", type=int, default=4)
    p_detect.set_defaults(func=detect.run)

    # ---------------------------------------------------------------
    # 2) build-controls
    # ---------------------------------------------------------------
    p_controls = subparsers.add_parser(
        "build-controls",
        help="Génère le JSON de référence (tailles max par outil) à partir des VCF fusionnés d'une cohorte de contrôles.",
    )
    p_controls.add_argument(
        "--controls-dir",
        required=True,
        help="Dossier de sortie de `detect` pour la cohorte contrôle "
        "({controls-dir}/{sample_id}/{sample_id}.merged.vcf).",
    )
    p_controls.add_argument(
        "--samples-list",
        help="Optionnel: TSV (colonne sample_id) pour restreindre aux échantillons listés "
        "(par défaut: tous les sous-dossiers de --controls-dir).",
    )
    p_controls.add_argument("-o", "--output", required=True, help="Chemin du JSON de sortie.")
    p_controls.set_defaults(func=controls.run)

    # ---------------------------------------------------------------
    # 3) compare
    # ---------------------------------------------------------------
    p_compare = subparsers.add_parser(
        "compare",
        help="Compare les expansions des patients au JSON de contrôles.",
    )
    p_compare.add_argument(
        "--patients-dir",
        required=True,
        help="Dossier de sortie de `detect` pour les patients "
        "({patients-dir}/{sample_id}/{sample_id}.merged.vcf).",
    )
    p_compare.add_argument(
        "--patients",
        nargs="+",
        help="Optionnel: liste d'identifiants patients à comparer "
        "(par défaut: tous les sous-dossiers de --patients-dir).",
    )
    p_compare.add_argument(
        "--controls-json", required=True, help="JSON généré par `build-controls`."
    )
    p_compare.add_argument(
        "--genes-bed", required=True, help="BED gzippé chrom/start/end/gene."
    )
    p_compare.add_argument(
        "--exons-bed", required=True, help="BED gzippé chrom/start/end/GENE_exonN (ex: MANE Select)."
    )
    p_compare.add_argument(
        "-x", "--threshold", type=int, default=0,
        help="Seuil additionnel: ne garder que diff_vs_controls > threshold (défaut: 0).",
    )
    p_compare.add_argument(
        "-t", "--triplet-only", action="store_true",
        help="Ne garder que les STR à motif de 3 pb ou plus.",
    )
    p_compare.add_argument(
        "-o", "--output", required=True, help="Fichier de sortie (TSV/CSV)."
    )
    p_compare.add_argument(
        "--format", choices=["tsv", "csv"], default="tsv", help="Format du rapport de sortie."
    )
    p_compare.set_defaults(func=compare.run)

    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args) or 0


if __name__ == "__main__":
    sys.exit(main())
