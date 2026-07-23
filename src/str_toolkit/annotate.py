"""
Annotation gène / feature (exon, UTR, intronique, intergénique, centromère,
télomère) pour un locus (chrom, pos). Porté de STRcompar2json.py.

Nécessite deux BED gzippés :
  - genes_bed : chrom, start, end, gene
  - exons_bed : chrom, start, end, "GENE_exonN"  (ex: exons MANE Select)
"""

from __future__ import annotations

import gzip
from pathlib import Path

CENTROMERE_COORDS = {
    "chr1": (121535434, 124535434),
    "chr2": (92326171, 95326171),
    "chr3": (90504854, 93504854),
    "chr4": (49660117, 52660117),
    "chr5": (46405641, 49405641),
    "chr6": (58626368, 61626368),
    "chr7": (58169654, 60828234),
    "chr8": (44033745, 45877265),
    "chr9": (43236168, 45518558),
    "chr10": (39686683, 41593521),
    "chr11": (51078349, 54425074),
    "chr12": (34769408, 37185252),
    "chr13": (16000001, 18051248),
    "chr14": (16000001, 18173523),
    "chr15": (17000001, 19725254),
    "chr16": (36311159, 38280682),
    "chr17": (22813680, 26885980),
    "chr18": (15460900, 20861206),
    "chr19": (24498981, 27190874),
    "chr20": (26436233, 30038348),
    "chr21": (10864561, 12915808),
    "chr22": (12954789, 15054318),
    "chrX": (58605580, 62412542),
    "chrY": (10316945, 10544039),
}

HG38_TELOMERES = {
    "chr1": {"p": (0, 10000), "q": (248946422, 248956422)},
    "chr2": {"p": (0, 10000), "q": (242183529, 242193529)},
    "chr3": {"p": (0, 10000), "q": (198285559, 198295559)},
    "chr4": {"p": (0, 10000), "q": (190204555, 190214555)},
    "chr5": {"p": (0, 10000), "q": (181528259, 181538259)},
    "chr6": {"p": (0, 10000), "q": (170795979, 170805979)},
    "chr7": {"p": (0, 10000), "q": (159335973, 159345973)},
    "chr8": {"p": (0, 10000), "q": (145128636, 145138636)},
    "chr9": {"p": (0, 10000), "q": (138384717, 138394717)},
    "chr10": {"p": (0, 10000), "q": (133787422, 133797422)},
    "chr11": {"p": (0, 10000), "q": (135076622, 135086622)},
    "chr12": {"p": (0, 10000), "q": (133265309, 133275309)},
    "chr13": {"p": (0, 10000), "q": (114354328, 114364328)},
    "chr14": {"p": (0, 10000), "q": (107033718, 107043718)},
    "chr15": {"p": (0, 10000), "q": (101981189, 101991189)},
    "chr16": {"p": (0, 10000), "q": (90328345, 90338345)},
    "chr17": {"p": (0, 10000), "q": (83247441, 83257441)},
    "chr18": {"p": (0, 10000), "q": (80363285, 80373285)},
    "chr19": {"p": (0, 10000), "q": (58607616, 58617616)},
    "chr20": {"p": (0, 10000), "q": (64434167, 64444167)},
    "chr21": {"p": (0, 10000), "q": (46699983, 46709983)},
    "chr22": {"p": (0, 10000), "q": (50808468, 50818468)},
    "chrX": {"p": (0, 10000), "q": None},
    "chrY": {"p": (0, 10000), "q": None},
}


def load_genes(path: str | Path) -> dict:
    dict_genes: dict = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            chrom, start, end, gene = line.strip().split("\t")
            dict_genes.setdefault(chrom, {})[gene] = {"start": int(start), "end": int(end)}
    return dict_genes


def load_exons(path: str | Path) -> dict:
    dict_exons: dict = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            chrom, start, end, gene_exon = line.strip().split("\t")
            gene, exon = gene_exon.split("_")
            dict_exons.setdefault(chrom, {}).setdefault(gene, {})[exon] = {
                "start": int(start),
                "end": int(end),
            }
    return dict_exons


def annotate_locus(chrom: str, pos: int, dict_genes: dict, dict_exons: dict) -> tuple[str, str]:
    """Retourne (genes_csv, features_csv) pour un locus donné, comme dans STRcompar2json.py."""
    list_genes: list[str] = []
    list_features: list[str] = []

    for gene, coords in dict_genes.get(chrom, {}).items():
        if not (coords["start"] < pos < coords["end"]):
            continue
        list_genes.append(gene)

        exons = dict_exons.get(chrom, {}).get(gene)
        if not exons:
            list_features.append("intronic")
            continue

        sorted_exons = sorted(exons.keys(), key=lambda x: int(x.replace("exon", "")))
        first_start = exons[sorted_exons[0]]["start"]
        last_start = exons[sorted_exons[-1]]["start"]
        strand = "+" if first_start < last_start else "-"

        if (strand == "+" and pos < exons[sorted_exons[0]]["start"]) or (
            strand == "-" and pos > exons[sorted_exons[0]]["end"]
        ):
            list_features.append("5'UTR")
            continue
        if (strand == "+" and pos > exons[sorted_exons[-1]]["end"]) or (
            strand == "-" and pos < exons[sorted_exons[-1]]["start"]
        ):
            list_features.append("3'UTR")
            continue

        for exon_name, exon_coords in exons.items():
            if exon_coords["start"] < pos < exon_coords["end"]:
                list_features.append(exon_name)
        if not any(f.startswith("exon") or f.endswith("UTR") for f in list_features):
            list_features.append("intronic")

    centro = CENTROMERE_COORDS.get(chrom)
    if centro and centro[0] < pos < centro[1]:
        list_features.append("centromere")

    telo = HG38_TELOMERES.get(chrom, {})
    if telo.get("q") and telo["q"][0] < pos < telo["q"][1]:
        list_features.append("telomere_q")
    if telo.get("p") and telo["p"][0] < pos < telo["p"][1]:
        list_features.append("telomere_p")

    if not list_genes:
        list_genes.append("intergenic")
    if not list_features:
        list_features.append(".")

    return ",".join(list_genes), ",".join(sorted(set(list_features)))
