import gzip

from str_toolkit import merge
from str_toolkit.compare import build_comparison_table
from str_toolkit.controls import collect_control_calls


def _write_merged_vcf(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    return merge.write_merged_vcf(records, path)


def _locus(chrom, pos, motif, sizes_by_source):
    return {
        "chrom": chrom,
        "pos": pos,
        "end": pos + 40,
        "motif": motif,
        "sources": sorted(sizes_by_source),
        "sizes_by_source": sizes_by_source,
    }


def test_collect_control_calls_builds_per_tool_max(tmp_path):
    controls_dir = tmp_path / "controls"

    _write_merged_vcf(
        controls_dir / "ctrl1" / "ctrl1.merged.vcf",
        [_locus("chr1", 1000, "AAAG", {"vamos_hap1": 20, "trgt_allele1": 18})],
    )
    _write_merged_vcf(
        controls_dir / "ctrl2" / "ctrl2.merged.vcf",
        [_locus("chr1", 1000, "AAAG", {"vamos_hap1": 25, "trgt_allele1": 15})],
    )

    registry = collect_control_calls(controls_dir)

    locus_id = "chr1_1000_AAAG"
    assert locus_id in registry
    assert registry[locus_id]["tools"]["vamos"]["max_size"] == 25
    assert registry[locus_id]["tools"]["vamos"]["n_observed"] == 2
    assert registry[locus_id]["tools"]["trgt"]["max_size"] == 18
    assert "tandem-genotypes" not in registry[locus_id]["tools"]


def _make_empty_bed_gz(path):
    with gzip.open(path, "wt") as fh:
        fh.write("")
    return path


def test_build_comparison_table_diff_per_tool_and_sort(tmp_path):
    patients_dir = tmp_path / "patients"

    # Locus 1 : expansion nette en VAMOS et TRGT
    # Locus 2 : légère expansion, un seul outil dispo -> diff plus petit
    _write_merged_vcf(
        patients_dir / "p01" / "p01.merged.vcf",
        [
            _locus("chr1", 1000, "AAAG", {"vamos_hap1": 60, "trgt_allele1": 55}),
            _locus("chr2", 5000, "AC", {"vamos_hap1": 12}),
        ],
    )

    controls_registry = {
        "chr1_1000_AAAG": {
            "chrom": "chr1", "pos": 1000, "motif": "AAAG",
            "tools": {"vamos": {"max_size": 20, "n_observed": 50}, "trgt": {"max_size": 18, "n_observed": 50}},
        },
        "chr2_5000_AC": {
            "chrom": "chr2", "pos": 5000, "motif": "AC",
            "tools": {"vamos": {"max_size": 10, "n_observed": 50}},
        },
    }

    genes_bed = _make_empty_bed_gz(tmp_path / "genes.bed.gz")
    exons_bed = _make_empty_bed_gz(tmp_path / "exons.bed.gz")

    df = build_comparison_table(
        patients_dir, ["p01"], controls_registry, str(genes_bed), str(exons_bed), threshold=0
    )

    assert len(df) == 2
    # Tri décroissant sur max_diff : locus1 (diff=40) avant locus2 (diff=2)
    assert df.iloc[0]["chrom"] == "chr1"
    assert df.iloc[0]["vamos_diff"] == 40
    assert df.iloc[0]["trgt_diff"] == 37
    assert df.iloc[0]["n_tools_expanded"] == 2
    assert df.iloc[0]["max_diff"] == 40

    assert df.iloc[1]["chrom"] == "chr2"
    assert df.iloc[1]["max_diff"] == 2
    assert df.iloc[1]["n_tools_expanded"] == 1
    assert df.iloc[1]["trgt_diff"] is None or df.iloc[1]["trgt_diff"] != df.iloc[1]["trgt_diff"]  # NaN/None


def test_build_comparison_table_filters_below_threshold(tmp_path):
    patients_dir = tmp_path / "patients"
    _write_merged_vcf(
        patients_dir / "p01" / "p01.merged.vcf",
        [_locus("chr1", 1000, "AAAG", {"vamos_hap1": 21})],
    )
    controls_registry = {
        "chr1_1000_AAAG": {
            "chrom": "chr1", "pos": 1000, "motif": "AAAG",
            "tools": {"vamos": {"max_size": 20, "n_observed": 50}},
        }
    }
    genes_bed = _make_empty_bed_gz(tmp_path / "genes.bed.gz")
    exons_bed = _make_empty_bed_gz(tmp_path / "exons.bed.gz")

    df = build_comparison_table(
        patients_dir, ["p01"], controls_registry, str(genes_bed), str(exons_bed), threshold=5
    )
    assert len(df) == 0
