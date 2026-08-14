import gzip

import pysam
import pytest

from str_toolkit import merge
from str_toolkit.instability import (
    _parse_allreads,
    compute_meiotic_instability,
    compute_somatic_instability,
    detect_mosaicism,
    is_sex_chromosome,
    match_transmitted_alleles,
    parse_longtr_for_somatic,
    read_duos,
    summarize_meiotic_instability_by_duo,
)


def _make_empty_bed_gz(path):
    with gzip.open(path, "wt") as fh:
        fh.write("")
    return path


def _write_merged_vcf(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    return merge.write_merged_vcf(records, path)


def _locus(chrom, pos, motif, sizes_by_source):
    return {
        "chrom": chrom, "pos": pos, "end": pos + 40, "motif": motif,
        "sources": sorted(sizes_by_source), "sizes_by_source": sizes_by_source,
    }


# ---------------------------------------------------------------------
# Meiotic instability
# ---------------------------------------------------------------------

def test_read_duos_rejects_invalid_duo_type(tmp_path):
    duos_tsv = tmp_path / "duos.tsv"
    duos_tsv.write_text("duo_id\tparent_id\tchild_id\tduo_type\nD1\tP1\tC1\tgrandmother_grandson\n")
    with pytest.raises(SystemExit):
        read_duos(str(duos_tsv))


def test_read_duos_accepts_valid_types(tmp_path):
    duos_tsv = tmp_path / "duos.tsv"
    duos_tsv.write_text("duo_id\tparent_id\tchild_id\tduo_type\nD1\tP1\tC1\tmother_son\n")
    duos = read_duos(str(duos_tsv))
    assert duos[0]["duo_type"] == "mother_son"


def test_match_transmitted_alleles_nearest_size():
    matches = match_transmitted_alleles(parent_sizes=[20.0, 40.0], child_sizes=[22.0, 41.0])
    matched_pairs = {(c, p) for c, p, _ in matches}
    assert (22.0, 20.0) in matched_pairs
    assert (41.0, 40.0) in matched_pairs
    diffs = {round(d, 1) for _, _, d in matches}
    assert diffs == {2.0, 1.0}


def test_match_transmitted_alleles_empty_inputs():
    assert match_transmitted_alleles([], [10.0]) == []
    assert match_transmitted_alleles([10.0], []) == []


def test_compute_meiotic_instability_end_to_end(tmp_path):
    data_dir = tmp_path / "data"
    _write_merged_vcf(
        data_dir / "parent1" / "parent1.merged.vcf",
        [_locus("chr1", 5000, "AAAG", {"vamos_hap1": 20, "vamos_hap2": 22})],
    )
    _write_merged_vcf(
        data_dir / "child1" / "child1.merged.vcf",
        [_locus("chr1", 5000, "AAAG", {"vamos_hap1": 21, "vamos_hap2": 30})],
    )

    duos = [{"duo_id": "D1", "parent_id": "parent1", "child_id": "child1", "duo_type": "mother_son"}]
    genes_bed = _make_empty_bed_gz(tmp_path / "genes.bed.gz")
    exons_bed = _make_empty_bed_gz(tmp_path / "exons.bed.gz")

    df = compute_meiotic_instability(data_dir, duos, str(genes_bed), str(exons_bed))

    assert len(df) == 2  # 2 child alleles matched
    assert df.iloc[0]["duo_type"] == "mother_son"
    assert df.iloc[0]["tool"] == "vamos"
    assert df.iloc[0]["location_category"] == "subtelomeric"
    assert df.iloc[0]["motif_category"] == "tetranucleotide"
    # child allele 21 -> nearest parent allele 20 (diff=1); child allele 30 -> nearest parent 22 (diff=8)
    diffs = sorted(df["diff"].tolist())
    assert diffs == [1.0, 8.0]


def test_is_sex_chromosome():
    assert is_sex_chromosome("chrX") is True
    assert is_sex_chromosome("chrY") is True
    assert is_sex_chromosome("X") is True
    assert is_sex_chromosome("Y") is True
    assert is_sex_chromosome("chr1") is False
    assert is_sex_chromosome("chr22") is False


def test_meiotic_instability_excludes_sex_chromosomes_by_default(tmp_path):
    data_dir = tmp_path / "data"
    _write_merged_vcf(
        data_dir / "father" / "father.merged.vcf",
        [
            _locus("chr1", 5000, "AAAG", {"vamos_hap1": 20}),
            _locus("chrX", 5000, "AAAG", {"vamos_hap1": 20}),
            _locus("chrY", 5000, "AAAG", {"vamos_hap1": 20}),
        ],
    )
    _write_merged_vcf(
        data_dir / "son" / "son.merged.vcf",
        [
            _locus("chr1", 5000, "AAAG", {"vamos_hap1": 25}),
            _locus("chrX", 5000, "AAAG", {"vamos_hap1": 25}),
            _locus("chrY", 5000, "AAAG", {"vamos_hap1": 25}),
        ],
    )

    duos = [{"duo_id": "D1", "parent_id": "father", "child_id": "son", "duo_type": "father_son"}]
    genes_bed = _make_empty_bed_gz(tmp_path / "genes.bed.gz")
    exons_bed = _make_empty_bed_gz(tmp_path / "exons.bed.gz")

    df = compute_meiotic_instability(data_dir, duos, str(genes_bed), str(exons_bed))
    assert set(df["chrom"]) == {"chr1"}  # X and Y dropped by default


def test_meiotic_instability_can_include_sex_chromosomes(tmp_path):
    data_dir = tmp_path / "data"
    _write_merged_vcf(
        data_dir / "father" / "father.merged.vcf",
        [
            _locus("chr1", 5000, "AAAG", {"vamos_hap1": 20}),
            _locus("chrX", 5000, "AAAG", {"vamos_hap1": 20}),
        ],
    )
    _write_merged_vcf(
        data_dir / "son" / "son.merged.vcf",
        [
            _locus("chr1", 5000, "AAAG", {"vamos_hap1": 25}),
            _locus("chrX", 5000, "AAAG", {"vamos_hap1": 25}),
        ],
    )

    duos = [{"duo_id": "D1", "parent_id": "father", "child_id": "son", "duo_type": "father_son"}]
    genes_bed = _make_empty_bed_gz(tmp_path / "genes.bed.gz")
    exons_bed = _make_empty_bed_gz(tmp_path / "exons.bed.gz")

    df = compute_meiotic_instability(
        data_dir, duos, str(genes_bed), str(exons_bed), exclude_sex_chromosomes=False
    )
    assert set(df["chrom"]) == {"chr1", "chrX"}


def test_compute_meiotic_instability_no_shared_loci_returns_empty(tmp_path):
    data_dir = tmp_path / "data"
    _write_merged_vcf(data_dir / "parent1" / "parent1.merged.vcf", [_locus("chr1", 5000, "AAAG", {"vamos_hap1": 20})])
    _write_merged_vcf(data_dir / "child1" / "child1.merged.vcf", [_locus("chr2", 9000, "AC", {"vamos_hap1": 10})])

    duos = [{"duo_id": "D1", "parent_id": "parent1", "child_id": "child1", "duo_type": "father_daughter"}]
    genes_bed = _make_empty_bed_gz(tmp_path / "genes.bed.gz")
    exons_bed = _make_empty_bed_gz(tmp_path / "exons.bed.gz")

    df = compute_meiotic_instability(data_dir, duos, str(genes_bed), str(exons_bed))
    assert len(df) == 0


def test_summarize_meiotic_instability_by_duo_avoids_pseudoreplication():
    import pandas as pd

    df = pd.DataFrame([
        {"duo_id": "D1", "duo_type": "mother_son", "tool": "vamos", "location_category": "exonic",
         "motif_category": "trinucleotide", "diff": 1.0},
        {"duo_id": "D1", "duo_type": "mother_son", "tool": "vamos", "location_category": "exonic",
         "motif_category": "trinucleotide", "diff": 3.0},
        {"duo_id": "D2", "duo_type": "mother_son", "tool": "vamos", "location_category": "exonic",
         "motif_category": "trinucleotide", "diff": 5.0},
    ])
    summary = summarize_meiotic_instability_by_duo(df)
    assert len(summary) == 2  # one row per duo, not per locus
    d1_row = summary[summary["duo_id"] == "D1"].iloc[0]
    assert d1_row["n_loci"] == 2
    assert d1_row["median_diff"] == 2.0


def test_summarize_meiotic_instability_by_duo_handles_empty_df():
    import pandas as pd

    df = pd.DataFrame(columns=["duo_id", "duo_type", "tool", "location_category", "motif_category", "diff"])
    summary = summarize_meiotic_instability_by_duo(df)
    assert len(summary) == 0


# ---------------------------------------------------------------------
# Somatic instability
# ---------------------------------------------------------------------

def test_parse_allreads_basic():
    pairs = _parse_allreads("-8|31;4|39")
    assert pairs == [(-8.0, 31), (4.0, 39)]


def test_parse_allreads_excludes_sentinel():
    pairs = _parse_allreads("-999|14;-8|31;4|39")
    assert (-999.0, 14) not in pairs
    assert len(pairs) == 2


def test_parse_allreads_none_returns_empty():
    assert _parse_allreads(None) == []


def test_detect_mosaicism_no_off_allele_reads():
    allreads = [(-8.0, 31), (4.0, 39)]
    result = detect_mosaicism("AAAG", called_alleles=[-8.0, 4.0], allreads=allreads, min_off_allele_reads=3)
    assert result["is_mosaic"] is False
    assert result["off_allele_reads"] == 0
    assert result["total_reads"] == 70


def test_detect_mosaicism_flags_sufficient_off_allele_reads():
    # Motif AAAG (len 4): a 3rd cluster at +12 is >=1 motif unit away from both called alleles
    allreads = [(-8.0, 31), (4.0, 39), (12.0, 5)]
    result = detect_mosaicism("AAAG", called_alleles=[-8.0, 4.0], allreads=allreads, min_off_allele_reads=3)
    assert result["is_mosaic"] is True
    assert result["off_allele_reads"] == 5
    assert result["mosaic_fraction"] == pytest.approx(5 / 75)


def test_detect_mosaicism_below_min_read_threshold_not_flagged():
    # Only 2 off-allele reads, below the default min_off_allele_reads=3
    allreads = [(-8.0, 31), (4.0, 39), (12.0, 2)]
    result = detect_mosaicism("AAAG", called_alleles=[-8.0, 4.0], allreads=allreads, min_off_allele_reads=3)
    assert result["is_mosaic"] is False


def test_detect_mosaicism_small_deviation_within_motif_unit_not_flagged():
    # A read 2bp off (less than the 4bp motif unit) should NOT count as off-allele
    allreads = [(-8.0, 31), (4.0, 39), (6.0, 10)]
    result = detect_mosaicism("AAAG", called_alleles=[-8.0, 4.0], allreads=allreads, min_off_allele_reads=3)
    assert result["is_mosaic"] is False
    assert result["off_allele_reads"] == 0


def _write_longtr_vcf_with_allreads(path):
    header = pysam.VariantHeader()
    header.add_line('##INFO=<ID=MOTIF,Number=1,Type=String,Description="motif">')
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="end">')
    header.add_line('##FORMAT=<ID=GB,Number=1,Type=String,Description="bp diff vs reference">')
    header.add_line('##FORMAT=<ID=ALLREADS,Number=1,Type=String,Description="per-read bp diffs">')
    header.add_line('##contig=<ID=chr1,length=248956422>')
    header.add_sample("sample01")
    with pysam.VariantFile(str(path), "w", header=header) as vf:
        record = vf.new_record(
            contig="chr1", start=999, stop=1040, alleles=("N", "<STR>"),
            info={"MOTIF": "AAAG", "END": 1040},
        )
        record.samples["sample01"]["GB"] = "-8|4"
        record.samples["sample01"]["ALLREADS"] = "-8|31;4|39;12|5"
        vf.write(record)


def test_parse_longtr_for_somatic(tmp_path):
    vcf_path = tmp_path / "sample.longtr.vcf.gz"
    uncompressed = tmp_path / "sample.longtr.vcf"
    _write_longtr_vcf_with_allreads(uncompressed)
    pysam.tabix_compress(str(uncompressed), str(vcf_path))

    records = parse_longtr_for_somatic(vcf_path)
    assert len(records) == 1
    assert records[0]["motif"] == "AAAG"
    assert records[0]["called_alleles"] == [-8.0, 4.0]
    assert (12.0, 5) in records[0]["allreads"]


def test_compute_somatic_instability_end_to_end(tmp_path):
    detect_dir = tmp_path / "detect"
    sample_dir = detect_dir / "s01"
    sample_dir.mkdir(parents=True)

    uncompressed = tmp_path / "s01.longtr.vcf"
    _write_longtr_vcf_with_allreads(uncompressed)
    pysam.tabix_compress(str(uncompressed), str(sample_dir / "s01.longtr.vcf.gz"))

    genes_bed = _make_empty_bed_gz(tmp_path / "genes.bed.gz")
    exons_bed = _make_empty_bed_gz(tmp_path / "exons.bed.gz")

    df = compute_somatic_instability(detect_dir, ["s01"], str(genes_bed), str(exons_bed))

    assert len(df) == 1
    assert df.iloc[0]["tool"] == "longtr"
    assert df.iloc[0]["is_mosaic"] == True
    assert df.iloc[0]["location_category"] == "subtelomeric"


def test_parse_longtr_for_somatic_c9orf72_real_data(tmp_path):
    """
    Regression test on a real LongTR line for C9orf72 (chr9:27573455,
    GCCCCG). GT="1|2", GB="-6|6721" (a 6bp contraction allele and a
    ~1120-unit pathogenic expansion allele), ALLREADS="-11|1;-7|1;-6|6;
    -4|2;6721|1" (11 reads total, matching INFO/DP=11). All the small
    stutter-range reads (-11..-4) are within one motif unit (6bp) of the
    -6 called allele, so this locus should NOT be flagged mosaic --
    validates that detect_mosaicism does not mistake normal stutter noise
    around a called allele for somatic mosaicism.
    """
    header = pysam.VariantHeader()
    header.add_line('##INFO=<ID=MOTIF,Number=1,Type=String,Description="motif">')
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="end">')
    header.add_line('##FORMAT=<ID=GB,Number=1,Type=String,Description="bp diff vs reference">')
    header.add_line('##FORMAT=<ID=ALLREADS,Number=1,Type=String,Description="per-read bp diffs">')
    header.add_line('##contig=<ID=chr9,length=138394717>')
    header.add_sample("patient01")

    uncompressed = tmp_path / "c9orf72.vcf"
    with pysam.VariantFile(str(uncompressed), "w", header=header) as vf:
        record = vf.new_record(
            contig="chr9", start=27573454, stop=27573581,
            alleles=("N", "<STR>"),
            info={"MOTIF": "GCCCCG", "END": 27573581},
        )
        record.samples["patient01"]["GB"] = "-6|6721"
        record.samples["patient01"]["ALLREADS"] = "-11|1;-7|1;-6|6;-4|2;6721|1"
        vf.write(record)
    vcf_path = tmp_path / "c9orf72.vcf.gz"
    pysam.tabix_compress(str(uncompressed), str(vcf_path))

    records = parse_longtr_for_somatic(vcf_path)
    assert len(records) == 1
    assert records[0]["called_alleles"] == [-6.0, 6721.0]
    assert sorted(records[0]["allreads"]) == [(-11.0, 1), (-7.0, 1), (-6.0, 6), (-4.0, 2), (6721.0, 1)]

    metrics = detect_mosaicism(
        "GCCCCG", records[0]["called_alleles"], records[0]["allreads"], min_off_allele_reads=3
    )
    assert metrics["total_reads"] == 11  # matches INFO/DP=11 on the real line
    assert metrics["is_mosaic"] is False
    assert metrics["off_allele_reads"] == 0
