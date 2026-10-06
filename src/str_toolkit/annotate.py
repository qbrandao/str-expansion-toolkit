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

Both layers go through `LocationIndex`, built once with
`build_location_index`. Build it at the start of a run and pass it to every
call: without it, each locus rescans every gene of its chromosome and
re-sorts each candidate gene's exon list, which measured 16 times slower on
a 20000-gene annotation.

Requires two gzipped BED files:
  - genes_bed: chrom, start, end, gene
  - exons_bed: chrom, start, end, "GENE_exonN"  (e.g. MANE Select exons)
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

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

# GRCh38 primary assembly chromosome lengths, used for the subtelomeric
# window below. Kept separate from HG38_TELOMERES, which holds the assembly's
# telomere GAP intervals: those 10 kb blocks are unsequenced N padding, so no
# repeat can ever be called inside them. Using them as the subtelomeric
# definition made the category structurally unreachable, which is what a run
# on a real genome showed: zero loci, with the first mappable locus of every
# chromosome sitting just past the 10 kb mark.
GRCH38_CHROM_LENGTHS = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
    "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
    "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
    "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
    "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}

# Distance from either chromosome end counted as subtelomeric. A scientific
# choice, not a constant of nature: the literature uses spans from about
# 100 kb to 500 kb. State the value used in any methods section, keep it the
# same across cases and controls, and set it to 0 to fall back to the strict
# telomere gap, which yields no loci.
DEFAULT_SUBTELOMERE_BP = 500_000

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
    # the q gaps were missing, so the q end of both was never annotated
    "chrX": {"p": (0, 10000), "q": (156030895, 156040895)},
    "chrY": {"p": (0, 10000), "q": (57217415, 57227415)},
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


def _bed_lines(path: str | Path):
    """
    Yields (chrom, start, end, name) for each usable BED line.

    Tolerant on purpose: a single header line or an extra column in a
    sixty-thousand-line annotation would otherwise abort a whole run at load
    time. Comment, track and browser lines are skipped, only the first four
    columns are read, and anything unparsable is counted and reported rather
    than raised.
    """
    skipped = 0
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                skipped += 1
                continue
            try:
                yield parts[0], int(parts[1]), int(parts[2]), parts[3]
            except ValueError:
                skipped += 1
    if skipped:
        logger.warning("%s: %d line(s) skipped as unparsable", path, skipped)


def load_genes(path: str | Path) -> dict:
    dict_genes: dict = {}
    for chrom, start, end, gene in _bed_lines(path):
        dict_genes.setdefault(chrom, {})[gene] = {"start": start, "end": end}
    return dict_genes


def load_exons(path: str | Path) -> dict:
    dict_exons: dict = {}
    malformed = 0
    for chrom, start, end, gene_exon in _bed_lines(path):
        # rsplit, so a gene name containing an underscore does not raise
        if "_" not in gene_exon:
            malformed += 1
            continue
        gene, exon = gene_exon.rsplit("_", 1)
        dict_exons.setdefault(chrom, {}).setdefault(gene, {})[exon] = {
            "start": start,
            "end": end,
        }
    if malformed:
        logger.warning(
            "%s: %d name(s) without a '_exonN' suffix were skipped. The fourth "
            "column must look like GENE_exon1.", path, malformed
        )
    return dict_exons


# ---------------------------------------------------------------------
# Position index
#
# Both annotation layers ask the same question millions of times: which
# genes contain this position. Scanning every gene of the chromosome per
# locus, and re-sorting each candidate gene's exon list while doing it,
# costs billions of Python operations on a two-million-locus call set.
#
# Genes are therefore bucketed into fixed-size bins once, and everything
# that does not depend on the queried position (strand, TSS, exon
# intervals) is precomputed per gene. A lookup then touches only the few
# genes sharing the position's bin.
# ---------------------------------------------------------------------

_BIN_SIZE = 100_000


class LocationIndex:
    """Prepared gene and exon intervals, built once and queried per locus."""

    def __init__(self, dict_genes: dict, dict_exons: dict):
        self._bins: dict[tuple[str, int], list[dict]] = {}
        for chrom, genes in dict_genes.items():
            for gene, coords in genes.items():
                entry = self._prepare(chrom, gene, coords, dict_exons.get(chrom, {}).get(gene))
                first_bin = entry["start"] // _BIN_SIZE
                last_bin = entry["end"] // _BIN_SIZE
                for b in range(first_bin, last_bin + 1):
                    self._bins.setdefault((chrom, b), []).append(entry)
        # a stable order, so an overlap between two genes does not resolve
        # differently from one run to the next
        for entries in self._bins.values():
            entries.sort(key=lambda e: (e["start"], e["gene"]))

    @staticmethod
    def _prepare(chrom: str, gene: str, coords: dict, exons: dict | None) -> dict:
        entry = {
            "gene": gene,
            "start": coords["start"],
            "end": coords["end"],
            "exon_intervals": [],
            "exon_names": [],
            "strand": None,
            "tss": None,
            "first_start": None,
            "first_end": None,
            "last_start": None,
            "last_end": None,
        }
        if not exons:
            return entry

        sorted_exons, strand = _gene_strand(exons)
        first, last = exons[sorted_exons[0]], exons[sorted_exons[-1]]
        entry.update(
            strand=strand,
            tss=first["start"] if strand == "+" else first["end"],
            first_start=first["start"],
            first_end=first["end"],
            last_start=last["start"],
            last_end=last["end"],
            exon_intervals=sorted(
                (exons[name]["start"], exons[name]["end"]) for name in sorted_exons
            ),
            exon_names=[(name, exons[name]["start"], exons[name]["end"]) for name in sorted_exons],
        )
        return entry

    def genes_at(self, chrom: str, pos: int) -> list[dict]:
        bucket = self._bins.get((chrom, pos // _BIN_SIZE))
        if not bucket:
            return []
        return [e for e in bucket if e["start"] < pos < e["end"]]


# Built lazily for the plain-dict entry points below. The dicts are kept
# alive by the cache, so their id() cannot be recycled under a live entry.
_INDEX_CACHE: dict[tuple[int, int], tuple[dict, dict, "LocationIndex"]] = {}


def build_location_index(dict_genes: dict, dict_exons: dict) -> LocationIndex:
    """
    Returns a LocationIndex, reusing a cached one for the same pair of dicts.

    Call this once per run and pass the result around when annotating many
    loci; the cache merely keeps the plain-dict entry points fast.
    """
    key = (id(dict_genes), id(dict_exons))
    hit = _INDEX_CACHE.get(key)
    if hit is not None:
        return hit[2]
    index = LocationIndex(dict_genes, dict_exons)
    _INDEX_CACHE[key] = (dict_genes, dict_exons, index)
    return index


def _exon_rank(name: str) -> tuple[int, str]:
    """
    Sort key for an exon name. 'exon12' sorts by its number; a name without a
    usable number falls back to lexical order after the numbered ones, rather
    than raising and taking the whole run down.
    """
    digits = "".join(c for c in name if c.isdigit())
    return (int(digits), "") if digits else (1 << 30, name)


def _gene_strand(exons: dict) -> tuple[list[str], str]:
    """Returns (exon names sorted 5'->3', strand) inferred from exon order."""
    sorted_exons = sorted(exons.keys(), key=_exon_rank)
    first_start = exons[sorted_exons[0]]["start"]
    last_start = exons[sorted_exons[-1]]["start"]
    strand = "+" if first_start < last_start else "-"
    return sorted_exons, strand


def annotate_locus(
    chrom: str,
    pos: int,
    dict_genes: dict,
    dict_exons: dict,
    index: LocationIndex | None = None,
) -> tuple[str, str]:
    """
    Returns (genes_csv, features_csv) for a given locus.

    Pass `index` (from build_location_index) when annotating many loci; it is
    otherwise derived from the two dicts and cached.
    """
    index = index or build_location_index(dict_genes, dict_exons)

    list_genes: list[str] = []
    list_features: list[str] = []

    for entry in index.genes_at(chrom, pos):
        list_genes.append(entry["gene"])

        if not entry["exon_names"]:
            list_features.append("intronic")
            continue

        strand = entry["strand"]
        if (strand == "+" and pos < entry["first_start"]) or (
            strand == "-" and pos > entry["first_end"]
        ):
            list_features.append("5'UTR")
            continue
        if (strand == "+" and pos > entry["last_end"]) or (
            strand == "-" and pos < entry["last_start"]
        ):
            list_features.append("3'UTR")
            continue

        hit_exon = False
        for exon_name, exon_start, exon_end in entry["exon_names"]:
            if exon_start < pos < exon_end:
                list_features.append(exon_name)
                hit_exon = True
        if not hit_exon:
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
    index: LocationIndex | None = None,
    subtelomere_bp: int = DEFAULT_SUBTELOMERE_BP,
) -> str:
    """
    Classifies a locus into ONE mutually exclusive genomic location category,
    for the genome-wide VNTR repertoire (paper Methods 2.8):

      - "subtelomeric"     : within subtelomere_bp of either chromosome end
      - "paracentromeric"  : within the centromeric window
      - "5prime_region"    : promoter window upstream of the TSS, or 5' UTR
      - "exonic"           : within an exon (3' UTR is counted as exonic,
                              since it is part of the terminal exon)
      - "intronic"         : within a gene but not exonic/5' region
      - "intergenic_other" : intergenic, and not subtelomeric/paracentromeric

    Priority: subtelomeric > paracentromeric > 5prime_region > exonic >
    intronic > intergenic_other. Where several genes overlap the position,
    the most specific category among them wins, so the result does not
    depend on the order the genes appear in the BED file.

    Pass `index` (from build_location_index) when classifying many loci.
    """
    # Measured from the chromosome ends rather than from the telomere gap,
    # see GRCH38_CHROM_LENGTHS.
    length = GRCH38_CHROM_LENGTHS.get(chrom)
    if length is not None and subtelomere_bp > 0:
        if pos <= subtelomere_bp or pos >= length - subtelomere_bp:
            return "subtelomeric"
    else:
        telo = HG38_TELOMERES.get(chrom, {})
        if telo.get("p") and telo["p"][0] < pos < telo["p"][1]:
            return "subtelomeric"
        if telo.get("q") and telo["q"][0] < pos < telo["q"][1]:
            return "subtelomeric"

    centro = CENTROMERE_COORDS.get(chrom)
    if centro and centro[0] < pos < centro[1]:
        return "paracentromeric"

    index = index or build_location_index(dict_genes, dict_exons)

    best: str | None = None
    for entry in index.genes_at(chrom, pos):
        category = _gene_category(entry, pos, promoter_bp)
        if category == "5prime_region":
            return category                      # highest gene-level priority
        if best is None or (best == "intronic" and category == "exonic"):
            best = category

    return best or "intergenic_other"


def _gene_category(entry: dict, pos: int, promoter_bp: int) -> str:
    """The location category contributed by one overlapping gene."""
    if not entry["exon_intervals"]:
        return "intronic"

    strand, tss = entry["strand"], entry["tss"]
    promoter_start = tss - promoter_bp if strand == "+" else tss
    promoter_end = tss if strand == "+" else tss + promoter_bp
    if promoter_start < pos < promoter_end:
        return "5prime_region"

    if (strand == "+" and pos < entry["first_start"]) or (
        strand == "-" and pos > entry["first_end"]
    ):
        return "5prime_region"  # annotated 5' UTR, outside the promoter window itself

    for exon_start, exon_end in entry["exon_intervals"]:
        if exon_start < pos < exon_end:
            return "exonic"
    # Includes 3' UTR: not modeled separately (see docstring), and any
    # position within the gene body that fell through the checks above.
    return "intronic"


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
