import pysam
import pytest

from str_toolkit.merge import _merge_overlapping_tg_rows, _read_tandem_genotypes_rows, _split_two_alleles, parse_longtr, parse_tandem_genotypes, parse_trgt, parse_vamos

TG_SAMPLE_LINES = [
    "chr1\t231799889\t231799934\tTCCCTTCCTCCCTTCC\t2.8\t.\t259,267,269,271,272,318,318,319\t.\n",
    "chr1\t231799881\t231799940\tTCCCTTCTTCCCTTCCTCCCTTCC\t2.4\t.\t173,178,179,181,181,212,212,213\t.\n",
    "chr1\t231799879\t231799940\tCTTCCCTT\t7.5\t.\t518,535,538,542,544,636,636,638\t.\n",
    "chr1\t231799833\t231799921\tTTCCC\t18.8\t.\t829,856,861,867,871,1018,1018,1022\t.\n",
]


def test_split_two_alleles_finds_the_biggest_gap():
    values = [259, 267, 269, 271, 272, 318, 318, 319]
    low, high = _split_two_alleles(values)
    assert low == 269  # median of [259,267,269,271,272]
    assert high == 318  # median of [318,318,319]


def test_split_two_alleles_single_value():
    assert _split_two_alleles([100]) == (100, 100)


def test_parse_tandem_genotypes_dedups_overlapping_motif_candidates(tmp_path):
    tsv_path = tmp_path / "sample.tandem_genotypes.tsv"
    tsv_path.write_text("".join(TG_SAMPLE_LINES))

    calls = parse_tandem_genotypes(tsv_path)

    # 4 overlapping candidate lines, all with 8 reads -> the last one
    # (motif TTCCC) wins (max() keeps the last one on ties).
    assert len(calls) == 2
    assert {c.source for c in calls} == {"tandem_genotypes_allele1", "tandem_genotypes_allele2"}
    kept_motif = calls[0].motif
    assert kept_motif in {"TCCCTTCCTCCCTTCC", "TCCCTTCTTCCCTTCCTCCCTTCC", "CTTCCCTT", "TTCCC"}
    for c in calls:
        assert c.motif == kept_motif
        assert c.raw["n_reads"] == 8


def test_merge_overlapping_tg_rows_keeps_most_reads():
    rows = [
        {"chrom": "chr1", "start": 100, "end": 150, "motif": "AC", "values": [10]},
        {"chrom": "chr1", "start": 110, "end": 160, "motif": "ACAC", "values": [20, 21, 22]},
        {"chrom": "chr1", "start": 300, "end": 350, "motif": "AC", "values": [15, 16]},
    ]
    kept = _merge_overlapping_tg_rows(rows)
    assert len(kept) == 2  # (100-150 + 110-160 overlap) and (300-350 isolated)
    n_reads = sorted(len(r["values"]) for r in kept)
    assert n_reads == [2, 3]


def _write_trgt_vcf(path):
    header = pysam.VariantHeader()
    header.add_line('##INFO=<ID=MOTIFS,Number=.,Type=String,Description="motifs">')
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="end">')
    header.add_line('##FORMAT=<ID=AL,Number=.,Type=Integer,Description="allele lengths">')
    header.add_line('##contig=<ID=chr1,length=248956422>')
    header.add_sample("sample01")
    with pysam.VariantFile(str(path), "w", header=header) as vf:
        record = vf.new_record(
            contig="chr1", start=999, stop=1040,
            alleles=("N", "<STR>"),
            info={"MOTIFS": "AAAG", "END": 1040},
        )
        record.samples["sample01"]["AL"] = (40, 55)
        vf.write(record)


def test_parse_trgt_reads_motifs_and_al(tmp_path):
    vcf_path = tmp_path / "sample.trgt.vcf.gz"
    uncompressed = tmp_path / "sample.trgt.vcf"
    _write_trgt_vcf(uncompressed)
    pysam.tabix_compress(str(uncompressed), str(vcf_path))

    calls = parse_trgt(vcf_path)
    assert len(calls) == 2
    assert {c.source for c in calls} == {"trgt_allele1", "trgt_allele2"}
    assert all(c.motif == "AAAG" for c in calls)
    assert sorted(c.size for c in calls) == [40.0, 55.0]


def _write_longtr_vcf(path):
    header = pysam.VariantHeader()
    header.add_line('##INFO=<ID=MOTIF,Number=1,Type=String,Description="motif">')
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="end">')
    # LongTR packs per-allele values into a single '|'-joined string (Number=1,
    # Type=String), NOT a comma-separated numeric array -- confirmed on a real
    # LongTR output line (C9orf72 locus, GB="-6|6721").
    header.add_line('##FORMAT=<ID=GB,Number=1,Type=String,Description="bp diff vs reference">')
    header.add_line('##contig=<ID=chr1,length=248956422>')
    header.add_sample("sample01")
    with pysam.VariantFile(str(path), "w", header=header) as vf:
        record = vf.new_record(
            contig="chr1", start=999, stop=1040,
            alleles=("N", "<STR>"),
            info={"MOTIF": "AAAG", "END": 1040},
        )
        record.samples["sample01"]["GB"] = "0|12"
        vf.write(record)


def test_parse_longtr_reads_motif_and_gb(tmp_path):
    vcf_path = tmp_path / "sample.longtr.vcf.gz"
    uncompressed = tmp_path / "sample.longtr.vcf"
    _write_longtr_vcf(uncompressed)
    pysam.tabix_compress(str(uncompressed), str(vcf_path))

    calls = parse_longtr(vcf_path)
    assert len(calls) == 2
    assert {c.source for c in calls} == {"longtr_allele1", "longtr_allele2"}
    assert all(c.motif == "AAAG" for c in calls)
    assert sorted(c.size for c in calls) == [0.0, 12.0]


def test_parse_longtr_missing_file_returns_empty(tmp_path):
    assert parse_longtr(tmp_path / "does_not_exist.vcf.gz") == []


def test_parse_longtr_c9orf72_real_data(tmp_path):
    """
    Regression test built from a real LongTR output line (C9orf72 locus,
    chr9:27573455, GCCCCG hexanucleotide repeat -- the ALS/FTD pathogenic
    expansion locus). GB packs alleles as a single '|'-joined string
    ("-6|6721"), NOT a comma-separated numeric array as the VCF spec would
    normally imply for a multi-valued FORMAT field. This locks in the fix
    for a bug where naive iteration over GB silently iterated over string
    characters instead of the two allele values.
    """
    header = pysam.VariantHeader()
    header.add_line('##INFO=<ID=MOTIF,Number=1,Type=String,Description="motif">')
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="end">')
    header.add_line('##FORMAT=<ID=GB,Number=1,Type=String,Description="bp diff vs reference">')
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
        vf.write(record)
    vcf_path = tmp_path / "c9orf72.vcf.gz"
    pysam.tabix_compress(str(uncompressed), str(vcf_path))

    calls = parse_longtr(vcf_path)
    assert len(calls) == 2
    assert all(c.motif == "GCCCCG" for c in calls)
    sizes = sorted(c.size for c in calls)
    assert sizes == [-6.0, 6721.0]  # a 6bp contraction allele and a ~1120-unit pathogenic expansion


def test_parse_longtr_second_real_locus_confirms_pipe_format(tmp_path):
    """
    Second independent real LongTR line (different locus, chr3:194943234,
    28bp motif, all-positive GB="1004|1350", 21 reads split 10/11 between
    strands). Confirms the pipe-separated GB/ALLREADS format is consistent
    across loci, not a one-off for the first (C9orf72) regression test.
    """
    header = pysam.VariantHeader()
    header.add_line('##INFO=<ID=MOTIF,Number=1,Type=String,Description="motif">')
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="end">')
    header.add_line('##FORMAT=<ID=GB,Number=1,Type=String,Description="bp diff vs reference">')
    header.add_line('##contig=<ID=chr3,length=198295559>')
    header.add_sample("patient01")

    uncompressed = tmp_path / "locus2.vcf"
    with pysam.VariantFile(str(uncompressed), "w", header=header) as vf:
        record = vf.new_record(
            contig="chr3", start=194943233, stop=194943717,
            alleles=("N", "<STR>"),
            info={"MOTIF": "CCACACTCTCCCACACTCTCCCACTCTC", "END": 194943717},
        )
        record.samples["patient01"]["GB"] = "1004|1350"
        vf.write(record)
    vcf_path = tmp_path / "locus2.vcf.gz"
    pysam.tabix_compress(str(uncompressed), str(vcf_path))

    calls = parse_longtr(vcf_path)
    assert len(calls) == 2
    assert sorted(c.size for c in calls) == [1004.0, 1350.0]


def _write_vamos_c9orf72_vcf(path):
    header = pysam.VariantHeader()
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="end">')
    header.add_line('##INFO=<ID=RU,Number=1,Type=String,Description="repeat units">')
    header.add_line('##INFO=<ID=LEN_H1,Number=1,Type=Integer,Description="length in motif units">')
    header.add_line('##contig=<ID=chr9,length=138394717>')
    header.add_sample("patient01")
    with pysam.VariantFile(str(path), "w", header=header) as vf:
        record = vf.new_record(
            contig="chr9", start=27573414, stop=27573546,
            alleles=("N", "<VNTR>"),
            info={
                "END": 27573546,
                "RU": "GGCCCC,GCCCC,GGGCCC,GCCCCC,GGCACCGC,AACCGC,TCACTC,ACCCACTC,"
                      "GCCACC,TGCGCC,GCGCCTCC,GCGCGCC,GGCGCA,GACCAC",
                "LEN_H1": 1123,
            },
        )
        vf.write(record)


def test_parse_vamos_c9orf72_real_data(tmp_path):
    """
    Regression test from a real VAMOS line for C9orf72 (chr9:27573415).
    RU lists 14 candidate motifs (comma-separated) -- only the first
    ("GGCCCC") is used as the locus motif. LEN_H1=1123 (repeat-motif
    units, not bp) is consistent with the LongTR regression test for the
    same gene (GB=6721bp / 6bp motif ~= 1120 units) -- cross-tool sanity
    check that VAMOS's unit convention (motif-repeat units) is correctly
    understood.
    """
    hap1_vcf = tmp_path / "patient01_assembly.hap1.vcf"
    _write_vamos_c9orf72_vcf(hap1_vcf)

    calls = parse_vamos({"hap1": hap1_vcf})
    assert len(calls) == 1
    assert calls[0].motif == "GGCCCC"
    assert calls[0].size == 1123.0
    assert calls[0].chrom == "chr9"


def test_parse_tandem_genotypes_pools_forward_and_reverse_strand_reads(tmp_path):
    """
    Regression test from a real tandem-genotypes line for C9orf72
    (chr9:27573484). Unlike earlier examples where column 8 (reverse
    strand) was always '.', this line has real data in BOTH column 7
    (forward strand) and column 8 (reverse strand) -- confirms both must
    be parsed and pooled, not just column 7 alone.
    """
    tsv_path = tmp_path / "sample.tandem_genotypes.tsv"
    tsv_path.write_text(
        "chr9\t27573484\t27573546\tGCCCCG\t10.8\t.\t"
        "-1,-1,-1,-1,-1,0,2,2,2,2,2\t-2,-1,-1,-1,2,2,2,2,2,3\n"
    )

    rows = _read_tandem_genotypes_rows(tsv_path)
    assert len(rows) == 1
    assert len(rows[0]["values"]) == 21  # 11 forward + 10 reverse reads pooled
