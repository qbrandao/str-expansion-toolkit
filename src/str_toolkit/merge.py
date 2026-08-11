"""
Merges VAMOS / TRGT / tandem-genotypes / LongTR outputs into a single VCF
per sample.

PROBLEM ADDRESSED: each tool anchors its coordinates differently for what
is biologically the same STR locus (VAMOS = assembly coordinates, TRGT =
catalog BED coordinates, tandem-genotypes = TRF coordinates, LongTR =
its own regions BED). Merging by strict equality of (chrom, pos, motif) is
therefore not possible: a tolerant interval + canonical motif matching
approach is required.

IMPORTANT on units: sizes are NOT necessarily comparable across tools
(VAMOS: length in motif-repeat units; TRGT and tandem-genotypes: length in
bp, but measured differently -- TRGT via direct genotyping, tandem-genotypes
via read-level length clustering; LongTR: bp difference from the reference,
not an absolute length). This module therefore does NOT merge sizes into a
single value: each merged locus keeps one size per source
(SIZES=vamos_hap1:42,trgt_allele1:38,...). Interpretation/comparison to
controls happens downstream, source by source (see str_toolkit/compare.py).
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median

import pysam

# ---------------------------------------------------------------------
# Normalized representation of an STR call, independent of the source tool
# ---------------------------------------------------------------------

@dataclass
class STRCall:
    chrom: str
    start: int
    end: int
    motif: str
    size: float
    source: str  # e.g. "vamos_hap1", "trgt_allele2", "tandem_genotypes_allele1"
    raw: dict = field(default_factory=dict)


def canonical_motif(motif: str) -> str:
    """Minimal circular rotation: AAGA / AGAA / GAAA / AAAG -> AAAG."""
    m = motif.strip().upper()
    if not m:
        return m
    rotations = [m[i:] + m[:i] for i in range(len(m))]
    return min(rotations)


# ---------------------------------------------------------------------
# Per-tool parsers -- adjust field names if your VCF/TSV files differ
# ---------------------------------------------------------------------

def _first(value):
    """INFO can be a single value or a tuple depending on Number= in the header."""
    if isinstance(value, (tuple, list)):
        return value[0] if value else None
    return value


def parse_vamos(hap_vcfs: dict[str, Path]) -> list[STRCall]:
    """
    Reads the per-haplotype VAMOS VCFs (e.g. {"hap1": ..., "hap2": ...}).
    Fields used (same as STRlist2json.py): INFO/RU (motif, first element if
    a list), INFO/LEN_H1 (size, in motif-repeat units).
    """
    calls = []
    for hap, vcf_path in hap_vcfs.items():
        vcf_path = Path(vcf_path)
        if not vcf_path.exists():
            continue
        with pysam.VariantFile(str(vcf_path)) as vf:
            for record in vf:
                ru = record.info.get("RU")
                length = record.info.get("LEN_H1")
                if ru is None or length is None:
                    continue
                motif = str(_first(ru)).split(",")[0]
                end = record.info.get("END", record.stop)
                calls.append(
                    STRCall(
                        chrom=record.chrom,
                        start=record.pos,
                        end=int(end) if end else record.pos,
                        motif=motif,
                        size=float(_first(length)),
                        source=f"vamos_{hap}",
                        raw={"info": dict(record.info)},
                    )
                )
    return calls


def parse_trgt(vcf_path: Path) -> list[STRCall]:
    """
    Reads the VCF produced by `trgt genotype` (.vcf.gz).

    NOTE: the field names below (INFO/MOTIFS or INFO/RU, FORMAT/AL) are the
    ones generally used by TRGT, but may vary by version. Check with
    `zcat sample.trgt.vcf.gz | grep '^##'` and adjust if needed.
    """
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        return []
    calls = []
    with pysam.VariantFile(str(vcf_path)) as vf:
        for record in vf:
            motif_field = record.info.get("MOTIFS") or record.info.get("RU")
            if motif_field is None:
                continue
            motif = str(_first(motif_field)).split(",")[0]
            end = record.info.get("END", record.stop)

            for sample_name, sample_data in record.samples.items():
                al = sample_data.get("AL")  # allele lengths in bp, e.g. (10, 14)
                if al is None:
                    continue
                for i, allele_len in enumerate(al):
                    if allele_len is None:
                        continue
                    calls.append(
                        STRCall(
                            chrom=record.chrom,
                            start=record.pos,
                            end=int(end) if end else record.pos,
                            motif=motif,
                            size=float(allele_len),
                            source=f"trgt_allele{i + 1}",
                            raw={"sample": sample_name},
                        )
                    )
    return calls


def _parse_pipe_values(raw) -> list[float]:
    """
    Parses a LongTR FORMAT value that packs multiple per-allele numbers into
    a single string joined by '|' (e.g. GB="-6|6721" for a heterozygous
    call), rather than the usual VCF comma-separated convention -- CONFIRMED
    on a real LongTR VCF line (C9orf72 locus, chr9:27573455). pysam returns
    this as a plain string (or a 1-tuple wrapping one) since the field is
    declared Number=1/Type=String in the header, not a numeric array --
    naive iteration over it (e.g. `for x in sample_data["GB"]`) silently
    iterates over characters instead of alleles. Always go through this
    function to parse GB.
    """
    if raw is None:
        return []
    if isinstance(raw, (tuple, list)):
        raw = "|".join(str(x) for x in raw if x is not None)
    values = []
    for chunk in str(raw).split("|"):
        try:
            values.append(float(chunk))
        except ValueError:
            continue
    return values


def parse_longtr(vcf_path: Path) -> list[STRCall]:
    """
    Reads the VCF produced by `LongTR` (--tr-vcf, bgzipped).

    Fields used (confirmed in the official gymrek-lab/LongTR README, and
    verified against a real output line): INFO/MOTIF (locus motif),
    INFO/END, FORMAT/GB (bp difference of each allele from the reference --
    NOT an absolute length, unlike TRGT/AL; packed as "allele1|allele2",
    see _parse_pipe_values). As with tandem-genotypes, this value remains
    internally consistent (control max / patient diff comparison within
    the same tool) even though it is not directly comparable to the
    absolute lengths reported by other tools.
    """
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        return []
    calls = []
    with pysam.VariantFile(str(vcf_path)) as vf:
        for record in vf:
            motif_field = record.info.get("MOTIF")
            if motif_field is None:
                continue
            motif = str(_first(motif_field)).split(",")[0]
            end = record.info.get("END", record.stop)

            for sample_name, sample_data in record.samples.items():
                gb = _parse_pipe_values(sample_data.get("GB"))  # bp diff vs reference per allele
                if not gb:
                    continue
                for i, bp_diff in enumerate(gb):
                    if bp_diff is None:
                        continue
                    calls.append(
                        STRCall(
                            chrom=record.chrom,
                            start=record.pos,
                            end=int(end) if end else record.pos,
                            motif=motif,
                            size=float(bp_diff),
                            source=f"longtr_allele{i + 1}",
                            raw={"sample": sample_name},
                        )
                    )
    return calls


def _split_two_alleles(values: list[float]) -> tuple[float, float]:
    """
    Splits a list of per-read lengths (bp) into 2 groups (short/long allele)
    at the largest gap between sorted values, then returns the median of
    each group. Robust to read-level noise, unlike a simple min/max.
    """
    values = sorted(values)
    if len(values) == 1:
        return values[0], values[0]
    gap_idx = max(range(len(values) - 1), key=lambda i: values[i + 1] - values[i])
    low, high = values[: gap_idx + 1], values[gap_idx + 1 :]
    return median(low), median(high)


def _merge_overlapping_tg_rows(rows: list[dict]) -> list[dict]:
    """
    repeats.trf.bed often lists several overlapping candidate motifs
    (different periods) for the same TRF locus -- tandem-genotypes then
    reports one row per candidate. Rows whose intervals overlap (per
    chromosome) are grouped, keeping only the one covered by the most
    reads as the locus representative (an unambiguous criterion, unlike
    column 5 -- see parse_tandem_genotypes).
    """
    by_chrom: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_chrom[r["chrom"]].append(r)

    kept = []
    for chrom_rows in by_chrom.values():
        chrom_rows.sort(key=lambda r: r["start"])
        current_group = []
        current_end = None
        for r in chrom_rows:
            if current_group and r["start"] <= current_end:
                current_group.append(r)
                current_end = max(current_end, r["end"])
            else:
                if current_group:
                    kept.append(max(current_group, key=lambda g: len(g["values"])))
                current_group = [r]
                current_end = r["end"]
        if current_group:
            kept.append(max(current_group, key=lambda g: len(g["values"])))
    return kept


def _read_tandem_genotypes_rows(tsv_path: Path) -> list[dict]:
    """
    Shared row-reading logic for tandem-genotypes TSVs: parses the 8-column
    format and deduplicates overlapping candidate motifs (see
    _merge_overlapping_tg_rows). Returns one dict per locus with the RAW
    per-read length list still intact (`values`) -- used both by
    parse_tandem_genotypes (which reduces `values` to 2 allele sizes) and
    by str_toolkit.instability (which needs the raw per-read distribution
    for somatic mosaicism detection).
    """
    tsv_path = Path(tsv_path)
    if not tsv_path.exists():
        return []

    rows = []
    with open(tsv_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 7:
                continue
            try:
                chrom, start, end, motif = cols[0], int(cols[1]), int(cols[2]), cols[3]
                values = [float(x) for x in cols[6].split(",") if x]
            except ValueError:
                continue
            if not values:
                continue
            rows.append({"chrom": chrom, "start": start, "end": end, "motif": motif, "values": values})

    return _merge_overlapping_tg_rows(rows)


def parse_tandem_genotypes(tsv_path: Path) -> list[STRCall]:
    """
    Reads the TSV produced by `tandem-genotypes repeats.bed alignments.maf`.

    8 tab-separated columns; column 5 (index 4, "gene name" / "score"
    depending on the official docs -- ambiguous on our files, values like
    2.8/18.8 that do not look like a gene name) is IGNORED, not used:
      0 chrom, 1 start, 2 end, 3 motif, 4 (ignored), 5 '.',
      6 per-read lengths (bp, comma-separated), 7 '.'

    Column 6 is a list of repeat lengths measured per individual read --
    it is split into 2 groups (short/long allele, see _split_two_alleles),
    whose median is taken as the allele size in bp (comparable to
    TRGT/AL). Overlapping candidate motifs (from repeats.trf.bed) are
    deduplicated, keeping the one covered by the most reads (see
    _merge_overlapping_tg_rows).
    """
    calls = []
    for row in _read_tandem_genotypes_rows(tsv_path):
        allele1, allele2 = _split_two_alleles(row["values"])
        for i, size in enumerate((allele1, allele2)):
            calls.append(
                STRCall(
                    chrom=row["chrom"],
                    start=row["start"],
                    end=row["end"],
                    motif=row["motif"],
                    size=size,
                    source=f"tandem_genotypes_allele{i + 1}",
                    raw={"n_reads": len(row["values"])},
                )
            )
    return calls


# ---------------------------------------------------------------------
# Fuzzy matching by interval (+ tolerance window) and canonical motif
# ---------------------------------------------------------------------

def cluster_calls(calls: list[STRCall], window: int = 25) -> list[list[STRCall]]:
    """
    Groups calls into loci: two calls belong to the same cluster if their
    [start, end] intervals overlap (with a `window` bp tolerance on either
    side) AND their canonical motifs are identical.

    Implemented as a single position-sorted pass (sweep), with a list of
    "open" clusters (still reachable given the tolerance) -- no O(n^2)
    comparison across the whole chromosome.
    """
    by_chrom: dict[str, list[STRCall]] = defaultdict(list)
    for c in calls:
        by_chrom[c.chrom].append(c)

    all_clusters: list[list[STRCall]] = []

    for chrom_calls in by_chrom.values():
        chrom_calls.sort(key=lambda c: c.start)
        open_clusters: list[dict] = []  # [{"end": int, "motifs": set, "calls": [...]}]

        for call in chrom_calls:
            cm = canonical_motif(call.motif)

            still_open = []
            for oc in open_clusters:
                if oc["end"] + window >= call.start:
                    still_open.append(oc)
                else:
                    all_clusters.append(oc["calls"])
            open_clusters = still_open

            match = next((oc for oc in open_clusters if cm in oc["motifs"]), None)
            if match is not None:
                match["calls"].append(call)
                match["end"] = max(match["end"], call.end)
                match["motifs"].add(cm)
            else:
                open_clusters.append({"end": call.end, "motifs": {cm}, "calls": [call]})

        all_clusters.extend(oc["calls"] for oc in open_clusters)

    return all_clusters


def build_locus_record(cluster: list[STRCall]) -> dict:
    chrom = cluster[0].chrom
    start = min(c.start for c in cluster)
    end = max(c.end for c in cluster)

    motif_counts = Counter(canonical_motif(c.motif) for c in cluster)
    motif = motif_counts.most_common(1)[0][0]

    sizes_by_source = {c.source: c.size for c in cluster}

    return {
        "chrom": chrom,
        "pos": start,
        "end": end,
        "motif": motif,
        "sources": sorted(sizes_by_source),
        "sizes_by_source": sizes_by_source,
    }


# ---------------------------------------------------------------------
# Writing the merged VCF
# ---------------------------------------------------------------------

_VCF_HEADER = """##fileformat=VCFv4.2
##source=str-expansion-toolkit merge_to_vcf
##INFO=<ID=MOTIF,Number=1,Type=String,Description="Canonical motif of the locus">
##INFO=<ID=END,Number=1,Type=Integer,Description="End of the locus (max across sources)">
##INFO=<ID=SOURCES,Number=.,Type=String,Description="Tools/haplotypes/alleles that reported this locus">
##INFO=<ID=SIZES,Number=.,Type=String,Description="source:size pairs (units differ by tool, see README)">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def write_merged_vcf(records: list[dict], out_path: Path) -> Path:
    records = sorted(records, key=lambda r: (r["chrom"], r["pos"]))
    with open(out_path, "w") as fh:
        fh.write(_VCF_HEADER)
        for r in records:
            sizes_str = ",".join(f"{src}:{val:g}" for src, val in r["sizes_by_source"].items())
            sources_str = ",".join(r["sources"])
            info = f"MOTIF={r['motif']};END={r['end']};SOURCES={sources_str};SIZES={sizes_str}"
            fh.write(f"{r['chrom']}\t{r['pos']}\t.\tN\t.\t.\t.\t{info}\n")
    return out_path


# ---------------------------------------------------------------------
# Entry point called by detect.merge_to_vcf
# ---------------------------------------------------------------------

def merge_tool_outputs(sample_id: str, tool_outputs: dict[str, object], outdir: Path, window: int = 25) -> Path:
    calls: list[STRCall] = []

    if "vamos" in tool_outputs:
        calls += parse_vamos(tool_outputs["vamos"])
    if "trgt" in tool_outputs:
        calls += parse_trgt(tool_outputs["trgt"])
    if "tandem-genotypes" in tool_outputs:
        calls += parse_tandem_genotypes(tool_outputs["tandem-genotypes"])
    if "longtr" in tool_outputs:
        calls += parse_longtr(tool_outputs["longtr"])

    clusters = cluster_calls(calls, window=window)
    records = [build_locus_record(c) for c in clusters]

    out_path = Path(outdir) / f"{sample_id}.merged.vcf"
    return write_merged_vcf(records, out_path)


# ---------------------------------------------------------------------
# Reading the merged VCF (used by build-controls / compare)
# ---------------------------------------------------------------------

# Source (e.g. "vamos_hap1", "tandem_genotypes_allele2") -> tool mapping
_TOOL_FAMILY_NAMES = {"vamos": "vamos", "trgt": "trgt", "tandem_genotypes": "tandem-genotypes", "longtr": "longtr"}


def tool_family(source: str) -> str:
    """'vamos_hap1' -> 'vamos'; 'trgt_allele1' -> 'trgt'; 'tandem_genotypes_allele1' -> 'tandem-genotypes'."""
    prefix = source.rsplit("_", 1)[0]
    return _TOOL_FAMILY_NAMES.get(prefix, prefix)


def parse_merged_vcf(path: Path):
    """
    Reads back a VCF produced by write_merged_vcf. Yields one dict per
    locus: {chrom, pos, end, motif, sizes_by_source}.
    """
    path = Path(path)
    if not path.exists():
        return
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, _id, _ref, _alt, _qual, _filt, info = fields[:8]

            info_dict = {}
            for kv in info.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info_dict[k] = v

            sizes_by_source = {}
            for pair in info_dict.get("SIZES", "").split(","):
                if ":" not in pair:
                    continue
                src, val = pair.rsplit(":", 1)
                try:
                    sizes_by_source[src] = float(val)
                except ValueError:
                    continue

            yield {
                "chrom": chrom,
                "pos": int(pos),
                "end": int(info_dict.get("END", pos)),
                "motif": info_dict.get("MOTIF", ""),
                "sizes_by_source": sizes_by_source,
            }
