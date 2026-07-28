import pysam
import pytest

from str_toolkit.merge import _merge_overlapping_tg_rows, _split_two_alleles, parse_longtr, parse_tandem_genotypes, parse_trgt

TG_SAMPLE_LINES = [
    "chr1\t231799889\t231799934\tTCCCTTCCTCCCTTCC\t2.8\t.\t259,267,269,271,272,318,318,319\t.\n",
    "chr1\t231799881\t231799940\tTCCCTTCTTCCCTTCCTCCCTTCC\t2.4\t.\t173,178,179,181,181,212,212,213\t.\n",
    "chr1\t231799879\t231799940\tCTTCCCTT\t7.5\t.\t518,535,538,542,544,636,636,638\t.\n",
    "chr1\t231799833\t231799921\tTTCCC\t18.8\t.\t829,856,861,867,871,1018,1018,1022\t.\n",
]


def test_split_two_alleles_finds_the_biggest_gap():
    values = [259, 267, 269, 271, 272, 318, 318, 319]
    low, high = _split_two_alleles(values)
    assert low == 269  # médiane de [259,267,269,271,272]
    assert high == 318  # médiane de [318,318,319]


def test_split_two_alleles_single_value():
    assert _split_two_alleles([100]) == (100, 100)


def test_parse_tandem_genotypes_dedups_overlapping_motif_candidates(tmp_path):
    tsv_path = tmp_path / "sample.tandem_genotypes.tsv"
    tsv_path.write_text("".join(TG_SAMPLE_LINES))

    calls = parse_tandem_genotypes(tsv_path)

    # 4 lignes candidates chevauchantes, toutes avec 8 reads -> la dernière
    # (motif TTCCC) l'emporte (max() garde le dernier en cas d'égalité).
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
    assert len(kept) == 2  # (100-150 + 110-160 se chevauchent) et (300-350 isolé)
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
    header.add_line('##FORMAT=<ID=GB,Number=.,Type=Integer,Description="bp diff vs reference">')
    header.add_line('##contig=<ID=chr1,length=248956422>')
    header.add_sample("sample01")
    with pysam.VariantFile(str(path), "w", header=header) as vf:
        record = vf.new_record(
            contig="chr1", start=999, stop=1040,
            alleles=("N", "<STR>"),
            info={"MOTIF": "AAAG", "END": 1040},
        )
        record.samples["sample01"]["GB"] = (0, 12)
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
