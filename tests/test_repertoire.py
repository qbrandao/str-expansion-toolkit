import gzip

from str_toolkit import merge
from str_toolkit.repertoire import build_repertoire, summarize_repertoire


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


def _make_empty_bed_gz(path):
    with gzip.open(path, "wt") as fh:
        fh.write("")
    return path


def test_build_repertoire_classifies_each_locus(tmp_path):
    controls_dir = tmp_path / "controls"

    # Subtelomeric locus, tetranucleotide motif
    _write_merged_vcf(
        controls_dir / "ctrl1" / "ctrl1.merged.vcf",
        [
            _locus("chr1", 5000, "AAAG", {"vamos_hap1": 20}),
            _locus("chr1", 50000000, "AC", {"vamos_hap1": 10}),  # intergenic_other, dinucleotide
        ],
    )

    genes_bed = _make_empty_bed_gz(tmp_path / "genes.bed.gz")
    exons_bed = _make_empty_bed_gz(tmp_path / "exons.bed.gz")

    df = build_repertoire(controls_dir, str(genes_bed), str(exons_bed))

    assert len(df) == 2
    row1 = df[df["pos"] == 5000].iloc[0]
    assert row1["location_category"] == "subtelomeric"
    assert row1["motif_category"] == "tetranucleotide"
    assert row1["vamos_max_size"] == 20

    row2 = df[df["pos"] == 50000000].iloc[0]
    assert row2["location_category"] == "intergenic_other"
    assert row2["motif_category"] == "dinucleotide"


def test_summarize_repertoire_counts_by_category(tmp_path):
    controls_dir = tmp_path / "controls"
    _write_merged_vcf(
        controls_dir / "ctrl1" / "ctrl1.merged.vcf",
        [
            _locus("chr1", 5000, "AAAG", {"vamos_hap1": 20}),
            _locus("chr1", 6000, "AAAT", {"vamos_hap1": 15}),
            _locus("chr1", 50000000, "AC", {"vamos_hap1": 10}),
        ],
    )
    genes_bed = _make_empty_bed_gz(tmp_path / "genes.bed.gz")
    exons_bed = _make_empty_bed_gz(tmp_path / "exons.bed.gz")

    df = build_repertoire(controls_dir, str(genes_bed), str(exons_bed))
    summary = summarize_repertoire(df)

    subtelo_tetra = summary[
        (summary["location_category"] == "subtelomeric") & (summary["motif_category"] == "tetranucleotide")
    ]
    assert subtelo_tetra["n_loci"].iloc[0] == 2  # both chr1:5000 and chr1:6000

    other_di = summary[
        (summary["location_category"] == "intergenic_other") & (summary["motif_category"] == "dinucleotide")
    ]
    assert other_di["n_loci"].iloc[0] == 1
