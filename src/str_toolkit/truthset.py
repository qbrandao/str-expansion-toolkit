"""
`validate-truthset` subcommand.

Validates the meiotic instability calls against the Platinum Pedigree
tandem repeat truthset for CEPH 1463 (Kronenberg et al., Nat Methods 2025),
available at:
  s3://platinum-pedigree-data/variants/tr_truthset/GRCh38/ceph_1463_tandem_repeats.oa.vcf.gz

WHAT THIS VALIDATES, AND WHY

The truthset is not a list of de novo mutations. It is a set of tandem
repeat genotypes that were retained because they segregate consistently
with Mendelian inheritance across the pedigree. At those loci, transmission
is therefore known to be consistent, so a parent-to-child size difference
called by our pipeline is, with high probability, a false positive rather
than genuine instability.

This gives an empirical false-positive rate for the nearest-size
transmitted-allele assignment of instability.py, which is otherwise
assumed rather than measured (see paper Methods 2.10 and Limitations 4.2).

Loci are matched to the truthset with the same tolerance used elsewhere in
the pipeline (interval overlap within a window, see merge.cluster_calls),
because the truthset and our callers do not share a coordinate convention.

Use `--inspect` first: it reports the INFO and FORMAT fields actually
present in your copy of the truthset, so the field names used for matching
can be confirmed rather than assumed.
"""

from __future__ import annotations

import logging
from collections import Counter
from pathlib import Path

import pandas as pd
import pysam

from str_toolkit.merge import canonical_motif

logger = logging.getLogger(__name__)

# INFO keys that may carry the repeat motif, in order of preference.
MOTIF_KEYS = ("MOTIF", "MOTIFS", "RU", "TRF_MOTIF", "UNIT")


def inspect_truthset(vcf_path: Path, n_records: int = 200) -> dict:
    """
    Reports the structure of a truthset VCF: samples, INFO/FORMAT keys
    declared in the header, and which keys are actually populated in the
    first records. Run this before trusting any field-name assumption.
    """
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        raise SystemExit(f"File not found: {vcf_path}")

    with pysam.VariantFile(str(vcf_path)) as vf:
        header_info = sorted(vf.header.info.keys())
        header_format = sorted(vf.header.formats.keys())
        samples = list(vf.header.samples)

        info_seen, format_seen = Counter(), Counter()
        examples = []
        for i, rec in enumerate(vf):
            if i >= n_records:
                break
            for k in rec.info.keys():
                info_seen[k] += 1
            for _, sd in rec.samples.items():
                for k in sd.keys():
                    format_seen[k] += 1
            if len(examples) < 3:
                examples.append({
                    "chrom": rec.chrom, "pos": rec.pos, "stop": rec.stop,
                    "ref_len": len(rec.ref) if rec.ref else None,
                    "n_alts": len(rec.alts) if rec.alts else 0,
                    "info": {k: str(rec.info[k])[:60] for k in list(rec.info.keys())[:8]},
                })

    return {
        "samples": samples,
        "header_info_keys": header_info,
        "header_format_keys": header_format,
        "info_keys_populated": dict(info_seen),
        "format_keys_populated": dict(format_seen),
        "example_records": examples,
        "motif_key_found": next((k for k in MOTIF_KEYS if k in info_seen), None),
    }


def load_truthset_loci(vcf_path: Path) -> list[dict]:
    """
    Reads truthset loci as {chrom, start, end, motif}. The motif is taken
    from the first INFO key in MOTIF_KEYS that is present, and is left
    empty when none is (locus matching then falls back to position only).
    """
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        raise SystemExit(f"File not found: {vcf_path}")

    loci = []
    with pysam.VariantFile(str(vcf_path)) as vf:
        for rec in vf:
            motif = ""
            for key in MOTIF_KEYS:
                if key in rec.info:
                    val = rec.info[key]
                    if isinstance(val, (tuple, list)):
                        val = val[0] if val else None
                    if val is None:
                        # pysam returns None for a valueless INFO field, which
                        # str() would turn into the literal string "None".
                        continue
                    candidate = str(val).split(",")[0].strip()
                    # "." is the VCF missing-value placeholder. Treating it as a
                    # real motif would make every locus fail motif matching,
                    # silently emptying the overlap.
                    if candidate and candidate != ".":
                        motif = candidate
                        break
            end = rec.info.get("END", rec.stop)
            loci.append({
                "chrom": rec.chrom,
                "start": rec.pos,
                "end": int(end) if end else rec.pos,
                "motif": motif,
            })
    return loci


def _index_by_chrom(loci: list[dict]) -> dict[str, list[dict]]:
    idx: dict[str, list[dict]] = {}
    for locus in loci:
        idx.setdefault(locus["chrom"], []).append(locus)
    for chrom in idx:
        idx[chrom].sort(key=lambda x: x["start"])
    return idx


def annotate_with_truthset(
    instability_df: pd.DataFrame,
    truthset_loci: list[dict],
    window: int = 25,
    require_motif_match: bool = True,
) -> pd.DataFrame:
    """
    Adds an `in_truthset` column to a per-locus meiotic instability table,
    matching on interval overlap within `window` bp and, when the truthset
    carries motifs, on canonical motif identity.
    """
    if instability_df.empty:
        out = instability_df.copy()
        out["in_truthset"] = pd.Series(dtype=bool)
        return out

    idx = _index_by_chrom(truthset_loci)
    truth_has_motifs = any(locus["motif"] for locus in truthset_loci)
    if require_motif_match and not truth_has_motifs:
        logger.warning(
            "Truthset carries no motif field, falling back to position-only matching. "
            "Run --inspect to confirm which INFO keys are present."
        )

    flags = []
    for row in instability_df.itertuples(index=False):
        candidates = idx.get(row.chrom, [])
        pos = int(row.pos)
        motif_c = canonical_motif(str(row.motif))
        hit = False
        for locus in candidates:
            if locus["end"] + window < pos:
                continue
            if locus["start"] - window > pos:
                break
            if require_motif_match and truth_has_motifs and locus["motif"]:
                if canonical_motif(locus["motif"]) != motif_c:
                    continue
            hit = True
            break
        flags.append(hit)

    out = instability_df.copy()
    out["in_truthset"] = flags
    return out


def false_positive_summary(df: pd.DataFrame, tolerance: float = 0.0) -> pd.DataFrame:
    """
    At truthset loci, transmission is known to be Mendelian-consistent, so
    a nonzero parent-to-child difference is most likely a false positive.
    Returns, per tool, the fraction of truthset loci at which a difference
    greater than `tolerance` is called.

    This is an upper bound on the false-positive rate rather than an exact
    one: the truthset reflects consistency of the reported genotypes, and a
    locus may still be genuinely unstable while passing pedigree filters.
    """
    truth = df[df["in_truthset"]]
    if truth.empty:
        return pd.DataFrame(columns=["tool", "n_truthset_loci", "n_nonzero_diff", "apparent_fp_rate"])

    rows = []
    for tool, sub in truth.groupby("tool"):
        n = len(sub)
        n_nonzero = int((sub["diff"].abs() > tolerance).sum())
        rows.append({
            "tool": tool,
            "n_truthset_loci": n,
            "n_nonzero_diff": n_nonzero,
            "apparent_fp_rate": n_nonzero / n if n else float("nan"),
        })
    return pd.DataFrame(rows).sort_values("tool").reset_index(drop=True)


def run(args) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    if args.inspect:
        info = inspect_truthset(args.truthset)
        print("Samples in truthset:", info["samples"] or "(none, sites-only VCF)")
        print("\nINFO keys declared in header:", ", ".join(info["header_info_keys"]) or "(none)")
        print("FORMAT keys declared in header:", ", ".join(info["header_format_keys"]) or "(none)")
        print("\nINFO keys populated in first records:")
        for k, v in sorted(info["info_keys_populated"].items()):
            print(f"   {k}: {v}")
        print("\nFORMAT keys populated in first records:")
        for k, v in sorted(info["format_keys_populated"].items()):
            print(f"   {k}: {v}")
        print("\nMotif key used for matching:", info["motif_key_found"] or "NONE FOUND (position-only matching)")
        print("\nExample records:")
        for ex in info["example_records"]:
            print("  ", ex)
        return 0

    if not args.instability:
        raise SystemExit("--instability is required unless --inspect is used")

    sep = "," if args.format == "csv" else "\t"
    df = pd.read_csv(args.instability, sep=sep)
    required = {"chrom", "pos", "motif", "tool", "diff"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(
            f"{args.instability}: missing column(s) {sorted(missing)}. "
            "Expected the per-locus output of `str-toolkit meiotic-instability`."
        )

    truthset_loci = load_truthset_loci(args.truthset)
    logger.info("Truthset loci loaded: %d", len(truthset_loci))

    annotated = annotate_with_truthset(
        df, truthset_loci, window=args.window, require_motif_match=not args.ignore_motif
    )
    n_hit = int(annotated["in_truthset"].sum())
    logger.info(
        "Instability rows overlapping the truthset: %d of %d (%.1f%%)",
        n_hit, len(annotated), 100 * n_hit / len(annotated) if len(annotated) else 0,
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    annotated.to_csv(out_path, sep=sep, index=False)
    logger.info("Annotated table written: %s", out_path)

    summary = false_positive_summary(annotated, tolerance=args.tolerance)
    if args.summary:
        summary_path = Path(args.summary)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(summary_path, sep=sep, index=False)
        logger.info("False-positive summary written: %s", summary_path)

    print("\nApparent false-positive rate at truthset loci (per tool):")
    print(summary.to_string(index=False) if not summary.empty else "  (no overlapping loci)")
    return 0
