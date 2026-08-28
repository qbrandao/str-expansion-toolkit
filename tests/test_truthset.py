import pandas as pd
import pysam
import pytest

from str_toolkit.truthset import (
    annotate_with_truthset,
    false_positive_summary,
    inspect_truthset,
    load_truthset_loci,
)


def _write_truthset_vcf(path, records, motif_key="MOTIF"):
    header = pysam.VariantHeader()
    header.add_line(f'##INFO=<ID={motif_key},Number=1,Type=String,Description="motif">')
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="end">')
    header.add_line('##contig=<ID=chr1,length=248956422>')
    header.add_line('##contig=<ID=chr2,length=242193529>')
    with pysam.VariantFile(str(path), "w", header=header) as vf:
        for chrom, pos, end, motif in records:
            rec = vf.new_record(
                contig=chrom, start=pos - 1, stop=end,
                alleles=("N", "<STR>"),
                info={motif_key: motif, "END": end},
            )
            vf.write(rec)
    return path


def _instability_df(rows):
    return pd.DataFrame(rows, columns=["chrom", "pos", "motif", "tool", "diff"])


def test_load_truthset_loci(tmp_path):
    vcf = _write_truthset_vcf(tmp_path / "t.vcf", [("chr1", 1000, 1040, "AAAG")])
    loci = load_truthset_loci(vcf)
    assert len(loci) == 1
    assert loci[0]["chrom"] == "chr1"
    assert loci[0]["motif"] == "AAAG"


def test_load_truthset_loci_missing_file(tmp_path):
    with pytest.raises(SystemExit):
        load_truthset_loci(tmp_path / "nope.vcf")


def test_annotate_matches_within_window_and_motif(tmp_path):
    vcf = _write_truthset_vcf(tmp_path / "t.vcf", [("chr1", 1000, 1040, "AAAG")])
    loci = load_truthset_loci(vcf)

    df = _instability_df([
        ("chr1", 1010, "AAGA", "vamos", 0.0),   # shifted coords, rotated motif -> match
        ("chr1", 50000, "AAAG", "vamos", 2.0),  # far away -> no match
        ("chr1", 1010, "AC", "vamos", 1.0),     # right place, wrong motif -> no match
        ("chr2", 1010, "AAAG", "vamos", 0.0),   # wrong chromosome -> no match
    ])
    out = annotate_with_truthset(df, loci, window=25)
    assert out["in_truthset"].tolist() == [True, False, False, False]


def test_annotate_ignore_motif(tmp_path):
    vcf = _write_truthset_vcf(tmp_path / "t.vcf", [("chr1", 1000, 1040, "AAAG")])
    loci = load_truthset_loci(vcf)
    df = _instability_df([("chr1", 1010, "AC", "vamos", 1.0)])
    out = annotate_with_truthset(df, loci, window=25, require_motif_match=False)
    assert out["in_truthset"].tolist() == [True]


def test_annotate_empty_input_returns_column(tmp_path):
    vcf = _write_truthset_vcf(tmp_path / "t.vcf", [("chr1", 1000, 1040, "AAAG")])
    loci = load_truthset_loci(vcf)
    out = annotate_with_truthset(_instability_df([]), loci)
    assert "in_truthset" in out.columns
    assert len(out) == 0


def test_false_positive_summary_counts_nonzero_diffs():
    df = pd.DataFrame([
        # 4 truthset loci for vamos, 1 with a nonzero diff -> 25% apparent FP
        {"chrom": "chr1", "pos": 1, "motif": "A", "tool": "vamos", "diff": 0.0, "in_truthset": True},
        {"chrom": "chr1", "pos": 2, "motif": "A", "tool": "vamos", "diff": 0.0, "in_truthset": True},
        {"chrom": "chr1", "pos": 3, "motif": "A", "tool": "vamos", "diff": 0.0, "in_truthset": True},
        {"chrom": "chr1", "pos": 4, "motif": "A", "tool": "vamos", "diff": 5.0, "in_truthset": True},
        # non-truthset rows must not be counted
        {"chrom": "chr1", "pos": 5, "motif": "A", "tool": "vamos", "diff": 9.0, "in_truthset": False},
        # a second tool
        {"chrom": "chr1", "pos": 6, "motif": "A", "tool": "longtr", "diff": 3.0, "in_truthset": True},
    ])
    s = false_positive_summary(df)
    vamos = s[s["tool"] == "vamos"].iloc[0]
    assert vamos["n_truthset_loci"] == 4
    assert vamos["n_nonzero_diff"] == 1
    assert vamos["apparent_fp_rate"] == pytest.approx(0.25)
    longtr = s[s["tool"] == "longtr"].iloc[0]
    assert longtr["apparent_fp_rate"] == pytest.approx(1.0)


def test_false_positive_summary_tolerance():
    df = pd.DataFrame([
        {"chrom": "chr1", "pos": 1, "motif": "A", "tool": "vamos", "diff": 2.0, "in_truthset": True},
        {"chrom": "chr1", "pos": 2, "motif": "A", "tool": "vamos", "diff": 10.0, "in_truthset": True},
    ])
    s = false_positive_summary(df, tolerance=5.0)
    assert s.iloc[0]["n_nonzero_diff"] == 1  # only the 10.0 exceeds tolerance


def test_false_positive_summary_no_overlap_returns_empty():
    df = pd.DataFrame([
        {"chrom": "chr1", "pos": 1, "motif": "A", "tool": "vamos", "diff": 1.0, "in_truthset": False},
    ])
    assert false_positive_summary(df).empty


def test_inspect_reports_structure(tmp_path):
    vcf = _write_truthset_vcf(tmp_path / "t.vcf", [("chr1", 1000, 1040, "AAAG")])
    info = inspect_truthset(vcf)
    assert "MOTIF" in info["header_info_keys"]
    assert info["motif_key_found"] == "MOTIF"
    assert len(info["example_records"]) == 1


def test_load_truthset_falls_back_to_alternative_motif_key(tmp_path):
    # A truthset using RU instead of MOTIF must still be readable.
    vcf = _write_truthset_vcf(tmp_path / "t2.vcf", [("chr1", 1000, 1040, "AAAG")], motif_key="RU")
    loci = load_truthset_loci(vcf)
    assert loci[0]["motif"] == "AAAG"


def test_position_only_matching_when_truthset_has_no_motifs(tmp_path):
    vcf = _write_truthset_vcf(tmp_path / "t3.vcf", [("chr1", 1000, 1040, "")])
    loci = load_truthset_loci(vcf)
    df = _instability_df([("chr1", 1010, "AAAG", "vamos", 1.0)])
    out = annotate_with_truthset(df, loci, window=25, require_motif_match=True)
    assert out["in_truthset"].tolist() == [True]
