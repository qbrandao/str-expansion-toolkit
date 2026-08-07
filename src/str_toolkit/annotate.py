"""
Genomic annotation for STR/VNTR loci.

Two complementary annotation layers are provided:

1. `annotate_locus` -- gene name + feature (exon, UTR, intronic, intergenic,
   centromere, telomere) for a single locus. Used by `compare` for
   per-patient diagnostic reporting. Ported from the original
   STRcompar2json.py script.

2. `classify_location` / `classify_motif` -- mutually exclusive genomic
   location and motif-length categories used to build the genome-wide
   VNTR repertoire (stratified analyses of variability/instability by
   location and motif class).

Requires two gzipped BED files:
  - genes_bed: chrom, start, end, gene
  - exons_bed: chrom, start, end, "GENE_exonN"  (e.g. MANE Select exons)
"""

from __future__ import annotations

import gzip
from pathlib import Path

# Centromere and telomere coordinates (hg38). These currently define the
# immediate centromeric/telomeric boundaries reused from the original
# pipeline; widen them if a broader "para-centromeric"/"subtelomeric"
# definition is intended for the genome-wide repertoire (see paper Methods
# 2.8 -- this is a scientific choice to confirm, not just an engineering one).
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

# Default promoter window: distance (bp) upstream of the transcription
# start site considered part of the "5' region including promoter"
# location category. Adjust to match the convention cited in the paper
# (e.g. Ensembl regulatory build) once decided.
DEFAULT_PROMOTER_WINDOW_BP = 2000

MOTIF_LENGTH_CATEGORIES = {
    1: "mononucleotide",
    2: "dinucleotide",
    3: "trinucleotide",
    4: "tetranucleotide",
    5: "pentanucleotide",
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


def _gene_strand(exons: dict) -> tuple[list[str], str]:
    """Returns (exon names sorted 5'->3', strand) inferred from exon order."""
    sorted_exons = sorted(exons.keys(), key=lambda x: int(x.replace("exon", "")))
    first_start = exons[sorted_exons[0]]["start"]
    last_start = exons[sorted_exons[-1]]["start"]
    strand = "+" if first_start < last_start else "-"
    return sorted_exons, strand


def annotate_locus(chrom: str, pos: int, dict_genes: dict, dict_exons: dict) -> tuple[str, str]:
    """Returns (genes_csv, features_csv) for a given locus, as in STRcompar2json.py."""
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

        sorted_exons, strand = _gene_strand(exons)

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


def classify_location(
    chrom: str,
    pos: int,
    dict_genes: dict,
    dict_exons: dict,
    promoter_bp: int = DEFAULT_PROMOTER_WINDOW_BP,
) -> str:
    """
    Classifies a locus into ONE mutually exclusive genomic location category,
    for the genome-wide VNTR repertoire (paper Methods 2.8):

      - "subtelomeric"     : within the telomeric window of a chromosome end
      - "paracentromeric"  : within the centromeric window
      - "5prime_region"    : promoter window upstream of the TSS, or 5' UTR
      - "exonic"           : within an exon (3' UTR is counted as exonic,
                              since it is part of the terminal exon)
      - "intronic"         : within a gene but not exonic/5' region
      - "intergenic_other" : intergenic, and not subtelomeric/paracentromeric

    Checked in this priority order: subtelomeric > paracentromeric > gene
    overlap (5prime_region > exonic > intronic) > intergenic_other.
    """
    telo = HG38_TELOMERES.get(chrom, {})
    if telo.get("p") and telo["p"][0] < pos < telo["p"][1]:
        return "subtelomeric"
    if telo.get("q") and telo["q"][0] < pos < telo["q"][1]:
        return "subtelomeric"

    centro = CENTROMERE_COORDS.get(chrom)
    if centro and centro[0] < pos < centro[1]:
        return "paracentromeric"

    for gene, coords in dict_genes.get(chrom, {}).items():
        if not (coords["start"] < pos < coords["end"]):
            continue

        exons = dict_exons.get(chrom, {}).get(gene)
        if not exons:
            return "intronic"

        sorted_exons, strand = _gene_strand(exons)
        tss = exons[sorted_exons[0]]["start"] if strand == "+" else exons[sorted_exons[0]]["end"]
        promoter_start = tss - promoter_bp if strand == "+" else tss
        promoter_end = tss if strand == "+" else tss + promoter_bp

        if promoter_start < pos < promoter_end:
            return "5prime_region"
        if (strand == "+" and pos < exons[sorted_exons[0]]["start"]) or (
            strand == "-" and pos > exons[sorted_exons[0]]["end"]
        ):
            return "5prime_region"  # annotated 5' UTR, outside the promoter window itself

        for exon_coords in exons.values():
            if exon_coords["start"] < pos < exon_coords["end"]:
                return "exonic"
        # Includes 3' UTR: not modeled separately (see docstring), and any
        # position within the gene body that fell through the checks above.
        return "intronic"

    return "intergenic_other"


def classify_motif(motif: str) -> str:
    """
    Classifies a repeat motif into a length-based category (paper Methods 2.9):
    mononucleotide / dinucleotide / trinucleotide / tetranucleotide /
    pentanucleotide / hexanucleotide_or_longer.

    Expects a single motif string (compound/comma-separated motif lists are
    already reduced to their first element by each tool's parser in merge.py).
    """
    length = len(motif.strip())
    return MOTIF_LENGTH_CATEGORIES.get(length, "hexanucleotide_or_longer")
