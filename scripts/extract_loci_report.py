#!/usr/bin/env python3
"""
Turn a merged VCF into a short, readable report for a named set of genes.

Why this exists. `compare` answers "is this locus larger than in the
controls", and needs the control registry that `build-controls` produces. With
no control cohort processed yet, that question cannot be answered at all, and
the merged VCF on its own is two million rows with per-tool sizes in three
different units. Neither is something to hand a collaborator.

What this produces instead is the one thing that is defensible without
controls: at each locus overlapping a gene of interest, what each tool
measured, in its own unit and converted to motif units where that conversion
is defined, with the number of tools that saw the locus. It states plainly that
no locus is called expanded.

Units, which is the whole difficulty:
  VAMOS             length in motif-repeat units, per haplotype
  LongTR            base-pair difference from the reference allele, per allele
  tandem-genotypes  base-pair length from read-level clustering, per allele
  TRGT              absolute allele length in base pairs

Only the three base-pair units convert to motif units, by dividing by the motif
length. VAMOS is already in motif units. A reference repeat length would be
needed to turn a LongTR delta into an absolute size, and the merged VCF does
not carry one, so the delta is reported as a delta.

Usage:
  ./scripts/extract_loci_report.py \\
      --merged-vcf results/local/1312_DFT/1312_DFT.merged.vcf \\
      --genes-bed genes.bed.gz \\
      --genes HTT FMR1 C9orf72 ATXN1 ATXN2 ATXN3 DMPK \\
      --flank 5000 \\
      -o 1312_DFT_report.tsv

  # or every gene named in a file, one per line
  ./scripts/extract_loci_report.py ... --genes-file disease_genes.txt ...
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from str_toolkit.merge import parse_merged_vcf, tool_family  # noqa: E402

# Unit of the size each tool reports, and whether it converts to motif units.
TOOL_UNITS = {
    "vamos": ("motif units", "already in motif units"),
    "longtr": ("bp difference from the reference", "divide by motif length"),
    "tandem-genotypes": ("bp", "divide by motif length"),
    "trgt": ("bp, absolute allele length", "divide by motif length"),
}
IN_MOTIF_UNITS = {"vamos"}


def open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def load_gene_windows(genes_bed: str, wanted: set[str], flank: int) -> list[tuple]:
    """
    Returns [(chrom, start, end, gene)] for the wanted genes, padded by flank.

    Gene names are matched case insensitively, since annotation sources
    disagree on the case of names such as C9orf72 and C9ORF72.
    """
    wanted_lower = {g.lower() for g in wanted}
    found: dict[str, tuple] = {}
    with open_maybe_gzip(genes_bed) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end, name = parts[0], parts[1], parts[2], parts[3]
            key = name.lower()
            if key not in wanted_lower:
                continue
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                continue
            # a gene can appear on several lines; keep the widest span
            if key in found:
                _c, s0, e0, n0 = found[key]
                found[key] = (chrom, min(s0, start_i), max(e0, end_i), n0)
            else:
                found[key] = (chrom, start_i, end_i, name)

    windows = [
        (c, max(0, s - flank), e + flank, n) for c, s, e, n in found.values()
    ]
    missing = wanted_lower - set(found)
    if missing:
        print(f"WARNING: not found in {genes_bed}: {', '.join(sorted(missing))}",
              file=sys.stderr)
    return windows


def index_windows(windows: list[tuple]) -> dict[str, list[tuple]]:
    """
    Groups the windows by chromosome, so a merged VCF of two million rows is
    not scanned against every window in turn. A long gene list would otherwise
    dominate the runtime.
    """
    by_chrom: dict[str, list[tuple]] = {}
    for chrom, start, end, name in windows:
        by_chrom.setdefault(chrom, []).append((start, end, name))
    for entries in by_chrom.values():
        entries.sort()
    return by_chrom


def gene_at(chrom: str, pos: int, by_chrom: dict[str, list[tuple]]) -> str | None:
    for start, end, name in by_chrom.get(chrom, ()):
        if start > pos:
            break          # sorted by start, so no later window can match
        if pos <= end:
            return name
    return None


def build_rows(merged_vcf: str, windows: list[tuple]) -> list[dict]:
    by_chrom = index_windows(windows)
    rows = []
    for locus in parse_merged_vcf(Path(merged_vcf)):
        gene = gene_at(locus["chrom"], locus["pos"], by_chrom)
        if gene is None:
            continue

        motif = locus["motif"]
        motif_len = len(motif) if motif else 0
        sizes = locus["sizes_by_source"]

        per_tool: dict[str, list[float]] = {}
        for source, value in sizes.items():
            per_tool.setdefault(tool_family(source), []).append(value)

        row = {
            "gene": gene,
            "chrom": locus["chrom"],
            "pos": locus["pos"],
            "end": locus["end"],
            "motif": motif,
            "motif_length_bp": motif_len,
            "reference_span_bp": locus["end"] - locus["pos"],
            "n_tools": len(per_tool),
            "tools": ",".join(sorted(per_tool)),
        }

        for tool in ("vamos", "longtr", "tandem-genotypes", "trgt"):
            values = sorted(per_tool.get(tool, []), reverse=True)
            col = tool.replace("-", "_")
            row[f"{col}_raw"] = "|".join(f"{v:g}" for v in values) if values else ""
            if not values:
                row[f"{col}_motif_units"] = ""
            elif tool in IN_MOTIF_UNITS:
                row[f"{col}_motif_units"] = "|".join(f"{v:g}" for v in values)
            elif motif_len > 0:
                row[f"{col}_motif_units"] = "|".join(
                    f"{v / motif_len:.1f}" for v in values
                )
            else:
                row[f"{col}_motif_units"] = ""

        # Largest measurement expressed in motif units, for sorting only. It
        # mixes an absolute size with a delta, so it ranks candidates rather
        # than measuring anything.
        candidates = []
        for tool in ("vamos", "longtr", "tandem-genotypes", "trgt"):
            cell = row[f"{tool.replace('-', '_')}_motif_units"]
            candidates += [abs(float(x)) for x in cell.split("|") if x]
        row["largest_motif_units_any_tool"] = f"{max(candidates):.1f}" if candidates else ""

        rows.append(row)

    rows.sort(
        key=lambda r: (
            -float(r["largest_motif_units_any_tool"] or 0),
            r["gene"],
            r["chrom"],
            r["pos"],
        )
    )
    return rows


HEADER_NOTE = [
    "# Tandem repeat measurements at the requested genes, from a merged",
    "# single-sample call set. NO LOCUS IS CALLED EXPANDED: that requires a",
    "# control cohort, which is not available yet, so these are measurements",
    "# rather than findings.",
    "#",
    "# Size units differ by tool and are NOT comparable as raw numbers:",
    "#   vamos_raw             length in motif-repeat units, one value per haplotype",
    "#   longtr_raw            bp DIFFERENCE from the reference allele, per allele",
    "#   tandem_genotypes_raw  bp length from read-level clustering, per allele",
    "#   trgt_raw              absolute allele length in bp",
    "# The *_motif_units columns divide the bp values by the motif length so",
    "# they can be read side by side. A LongTR value stays a difference even",
    "# after conversion, so a positive number means longer than the reference.",
    "#",
    "# n_tools is how many of the tools reported this locus. Their catalogs",
    "# differ, so a locus seen by one tool only is common and is not evidence",
    "# of anything by itself.",
    "#",
    "# largest_motif_units_any_tool sorts the table. It mixes an absolute size",
    "# with a difference, so it ranks candidates and measures nothing.",
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--merged-vcf", required=True, help="{sample}.merged.vcf from detect.")
    ap.add_argument("--genes-bed", required=True, help="Gzipped or plain BED: chrom/start/end/gene.")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--genes", nargs="+", help="Gene names.")
    g.add_argument("--genes-file", help="File of gene names, one per line.")
    ap.add_argument("--flank", type=int, default=5000,
                    help="bp added on each side of every gene (default: 5000).")
    ap.add_argument("--min-tools", type=int, default=1,
                    help="Only keep loci reported by at least this many tools (default: 1).")
    ap.add_argument("-o", "--output", required=True, help="Output TSV.")
    ap.add_argument("--no-header-note", action="store_true",
                    help="Omit the explanatory comment block.")
    args = ap.parse_args()

    if args.genes_file:
        # one name per line, with an optional trailing comment after '#'
        wanted = set()
        for line in open(args.genes_file):
            name = line.split("#", 1)[0].strip()
            if name:
                wanted.add(name)
    else:
        wanted = set(args.genes)

    windows = load_gene_windows(args.genes_bed, wanted, args.flank)
    if not windows:
        print("ERROR: none of the requested genes were found in the BED.", file=sys.stderr)
        return 1
    print(f"{len(windows)} gene window(s), flank {args.flank} bp", file=sys.stderr)

    rows = [r for r in build_rows(args.merged_vcf, windows)
            if r["n_tools"] >= args.min_tools]
    if not rows:
        print("No locus of the merged VCF falls in those windows.", file=sys.stderr)

    fieldnames = [
        "gene", "chrom", "pos", "end", "motif", "motif_length_bp",
        "reference_span_bp", "n_tools", "tools",
        "largest_motif_units_any_tool",
        "vamos_raw", "vamos_motif_units",
        "longtr_raw", "longtr_motif_units",
        "tandem_genotypes_raw", "tandem_genotypes_motif_units",
        "trgt_raw", "trgt_motif_units",
    ]
    with open(args.output, "w", newline="") as fh:
        if not args.no_header_note:
            fh.write("\n".join(HEADER_NOTE) + "\n")
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"{len(rows)} locus row(s) -> {args.output}", file=sys.stderr)
    by_gene: dict[str, int] = {}
    for r in rows:
        by_gene[r["gene"]] = by_gene.get(r["gene"], 0) + 1
    for gene in sorted(by_gene):
        print(f"  {gene:<14} {by_gene[gene]:>5} loci", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
