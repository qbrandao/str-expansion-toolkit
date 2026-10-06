from str_toolkit.annotate import classify_location, classify_motif


def test_classify_motif_by_length():
    assert classify_motif("A") == "mononucleotide"
    assert classify_motif("AC") == "dinucleotide"
    assert classify_motif("AAG") == "trinucleotide"
    assert classify_motif("AAAG") == "tetranucleotide"
    assert classify_motif("AAAAG") == "pentanucleotide"
    assert classify_motif("AAAAAG") == "hexanucleotide_or_longer"
    assert classify_motif("AAAAAAAAAA") == "hexanucleotide_or_longer"


def test_classify_location_subtelomeric():
    assert classify_location("chr1", 5000, {}, {}) == "subtelomeric"
    assert classify_location("chr1", 248950000, {}, {}) == "subtelomeric"


def test_classify_location_paracentromeric():
    assert classify_location("chr1", 123000000, {}, {}) == "paracentromeric"


def test_classify_location_intergenic_other():
    assert classify_location("chr1", 50000000, {}, {}) == "intergenic_other"


def test_classify_location_exonic_and_intronic_and_5prime():
    genes = {"chr1": {"GENE1": {"start": 5001000, "end": 5002000}}}
    exons = {"chr1": {"GENE1": {"exon1": {"start": 5001100, "end": 5001200}, "exon2": {"start": 5001500, "end": 5001600}}}}

    assert classify_location("chr1", 5001150, genes, exons) == "exonic"
    assert classify_location("chr1", 5001300, genes, exons) == "intronic"
    # upstream of TSS (exon1 start = 5001100, plus strand since exon1 < exon2) -> promoter window
    assert classify_location("chr1", 5001050, genes, exons, promoter_bp=2000) == "5prime_region"


def test_classify_location_no_exons_defaults_to_intronic():
    genes = {"chr1": {"GENE1": {"start": 5001000, "end": 5002000}}}
    assert classify_location("chr1", 5001500, genes, {}) == "intronic"


# ---------------------------------------------------------------------
# The position index
# ---------------------------------------------------------------------

from str_toolkit.annotate import (  # noqa: E402
    LocationIndex,
    _gene_strand,
    annotate_locus,
    build_location_index,
    load_exons,
    load_genes,
)


def _two_genes():
    genes = {"chr1": {"GENEA": {"start": 5_000_000, "end": 5_010_000},
                      "GENEB": {"start": 5_005_000, "end": 5_020_000}}}
    exons = {"chr1": {
        # GENEA on the plus strand, two exons
        "GENEA": {"exon1": {"start": 5_001_000, "end": 5_001_200},
                  "exon2": {"start": 5_009_000, "end": 5_009_300}},
        # GENEB on the plus strand, its first exon starts inside GENEA
        "GENEB": {"exon1": {"start": 5_006_000, "end": 5_006_200},
                  "exon2": {"start": 5_019_000, "end": 5_019_400}},
    }}
    return genes, exons


def test_index_finds_every_overlapping_gene():
    genes, exons = _two_genes()
    index = build_location_index(genes, exons)

    names = sorted(e["gene"] for e in index.genes_at("chr1", 5_006_100))
    assert names == ["GENEA", "GENEB"]
    assert [e["gene"] for e in index.genes_at("chr1", 5_002_000)] == ["GENEA"]
    assert index.genes_at("chr1", 4_000_000) == []
    assert index.genes_at("chr2", 5_006_100) == []


def test_index_spans_genes_crossing_a_bin_boundary():
    # a gene wider than one bin must be found anywhere inside it
    genes = {"chr1": {"BIG": {"start": 1_000_000, "end": 1_450_000}}}
    index = build_location_index(genes, {})
    for pos in (1_000_001, 1_099_999, 1_100_001, 1_250_000, 1_449_999):
        assert [e["gene"] for e in index.genes_at("chr1", pos)] == ["BIG"]
    assert index.genes_at("chr1", 1_450_001) == []


def test_overlapping_genes_resolve_to_the_most_specific_category():
    """
    The documented priority is 5prime_region > exonic > intronic. With two
    overlapping genes the result must not depend on which one the BED lists
    first, which is what scanning until the first hit used to do.
    """
    genes, exons = _two_genes()
    # inside GENEA's body (intronic there) and inside GENEB's first exon
    assert classify_location("chr1", 5_006_100, genes, exons) == "exonic"
    # inside GENEA's body and in GENEB's promoter window
    assert classify_location("chr1", 5_005_100, genes, exons, promoter_bp=2000) == "5prime_region"

    # reversing the BED order must not change either answer
    reversed_genes = {"chr1": dict(reversed(list(genes["chr1"].items())))}
    assert classify_location("chr1", 5_006_100, reversed_genes, exons) == "exonic"
    assert classify_location("chr1", 5_005_100, reversed_genes, exons, promoter_bp=2000) == "5prime_region"


def test_classify_location_without_an_index_still_works():
    genes, exons = _two_genes()
    # no index argument: one is derived and cached
    assert classify_location("chr1", 5_002_000, genes, exons) == "intronic"
    assert classify_location("chr1", 5_001_100, genes, exons) == "exonic"


def test_build_location_index_is_cached_per_dict_pair():
    genes, exons = _two_genes()
    assert build_location_index(genes, exons) is build_location_index(genes, exons)
    other_genes, other_exons = _two_genes()
    assert build_location_index(other_genes, other_exons) is not build_location_index(genes, exons)


def test_annotate_locus_reports_all_overlapping_genes():
    genes, exons = _two_genes()
    gene_csv, feature_csv = annotate_locus("chr1", 5_006_100, genes, exons)
    assert gene_csv == "GENEA,GENEB"
    assert "exon1" in feature_csv      # GENEB exon1
    assert "intronic" in feature_csv   # GENEA


def test_annotate_locus_intergenic():
    genes, exons = _two_genes()
    assert annotate_locus("chr1", 4_000_000, genes, exons) == ("intergenic", ".")


def test_exon_name_with_an_underscore_in_the_gene(tmp_path):
    """load_exons must not choke on a gene name containing an underscore."""
    import gzip

    path = tmp_path / "exons.bed.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("chr1\t100\t200\tMT_ND1_exon1\n")
        fh.write("chr1\t300\t400\tMT_ND1_exon2\n")

    loaded = load_exons(path)
    assert set(loaded["chr1"]["MT_ND1"]) == {"exon1", "exon2"}


def test_minus_strand_promoter_is_upstream_in_genomic_terms():
    genes = {"chr1": {"MINUS": {"start": 2_000_000, "end": 2_050_000}}}
    # exon1 is the 5' end, so on the minus strand it has the HIGHER coordinate
    exons = {"chr1": {"MINUS": {
        "exon1": {"start": 2_049_000, "end": 2_049_500},
        "exon2": {"start": 2_001_000, "end": 2_001_500},
    }}}
    index = build_location_index(genes, exons)
    entry = index.genes_at("chr1", 2_025_000)[0]
    assert entry["strand"] == "-"
    assert entry["tss"] == 2_049_500
    # promoter sits above the TSS on the minus strand
    assert classify_location("chr1", 2_049_600, genes, exons, promoter_bp=2000, index=index) == "5prime_region"
    assert classify_location("chr1", 2_025_000, genes, exons, promoter_bp=2000, index=index) == "intronic"


# ---------------------------------------------------------------------
# BED loading, which must survive the small irregularities real annotation
# files carry rather than aborting a run at load time.
# ---------------------------------------------------------------------

def _write_gz(path, text):
    import gzip
    with gzip.open(path, "wt") as fh:
        fh.write(text)
    return path


def test_load_genes_tolerates_headers_extra_columns_and_junk(tmp_path, caplog):
    path = _write_gz(tmp_path / "genes.bed.gz",
                     "track name=genes\n"
                     "#chrom\tstart\tend\tgene\n"
                     "\n"
                     "chr1\t11121\t24894\tDDX11L16\n"
                     "chr1\t12010\t13670\tDDX11L1\t+\t0\n"      # extra columns
                     "chr1\tNOTANUMBER\t100\tBROKEN\n"           # unparsable
                     "chr1\t500\n")                              # too few columns

    with caplog.at_level("WARNING"):
        genes = load_genes(path)

    assert set(genes["chr1"]) == {"DDX11L16", "DDX11L1"}
    assert genes["chr1"]["DDX11L1"] == {"start": 12010, "end": 13670}
    assert "2 line(s) skipped" in caplog.text


def test_load_exons_reports_names_without_an_exon_suffix(tmp_path, caplog):
    path = _write_gz(tmp_path / "exons.bed.gz",
                     "chr1\t65418\t65433\tOR4F5_exon1\n"
                     "chr1\t65519\t65573\tOR4F5_exon2\n"
                     "chr1\t70000\t70100\tNOSUFFIX\n")

    with caplog.at_level("WARNING"):
        exons = load_exons(path)

    assert set(exons["chr1"]["OR4F5"]) == {"exon1", "exon2"}
    assert "NOSUFFIX" not in exons["chr1"]
    assert "1 name(s) without a '_exonN' suffix" in caplog.text


def test_exon_rank_orders_numerically_and_survives_odd_names():
    from str_toolkit.annotate import _exon_rank

    names = ["exon10", "exon2", "exon1", "exonA"]
    assert sorted(names, key=_exon_rank) == ["exon1", "exon2", "exon10", "exonA"]


def test_gene_strand_with_an_unnumbered_exon_does_not_raise():
    exons = {"exon1": {"start": 100, "end": 200},
             "exonX": {"start": 300, "end": 400}}
    names, strand = _gene_strand(exons)
    assert names[0] == "exon1"
    assert strand in ("+", "-")


def test_plain_uncompressed_bed_also_loads(tmp_path):
    path = tmp_path / "genes.bed"
    path.write_text("chr1\t100\t200\tGENEA\n")
    assert load_genes(path)["chr1"]["GENEA"] == {"start": 100, "end": 200}


# ---------------------------------------------------------------------
# The subtelomeric window
#
# The assembly's telomere annotation is a 10 kb gap of N at each chromosome
# end. No read aligns there, so using it as the subtelomeric definition made
# the category unreachable: a real genome produced zero subtelomeric loci,
# with the first callable locus of every chromosome sitting just past 10 kb.
# ---------------------------------------------------------------------

from str_toolkit.annotate import (  # noqa: E402
    DEFAULT_SUBTELOMERE_BP,
    GRCH38_CHROM_LENGTHS,
    HG38_TELOMERES,
)


def test_a_locus_just_past_the_telomere_gap_is_subtelomeric():
    # chr1's first callable locus in a real run was at 10757, 757 bp past the
    # 10 kb gap, and used to be classified intergenic_other
    assert classify_location("chr1", 10757, {}, {}) == "subtelomeric"


def test_both_chromosome_ends_are_covered():
    length = GRCH38_CHROM_LENGTHS["chr1"]
    assert classify_location("chr1", 1, {}, {}) == "subtelomeric"
    assert classify_location("chr1", length - 1, {}, {}) == "subtelomeric"
    assert classify_location("chr1", DEFAULT_SUBTELOMERE_BP + 1, {}, {}) != "subtelomeric"
    assert classify_location("chr1", length - DEFAULT_SUBTELOMERE_BP - 1, {}, {}) != "subtelomeric"


def test_sex_chromosome_q_ends_are_covered():
    """chrX and chrY carried q: None, so their q end was never annotated."""
    for chrom in ("chrX", "chrY"):
        length = GRCH38_CHROM_LENGTHS[chrom]
        assert HG38_TELOMERES[chrom]["q"] is not None
        assert classify_location(chrom, length - 1000, {}, {}) == "subtelomeric"


def test_window_is_configurable_and_zero_falls_back_to_the_gap():
    # a 100 kb window excludes a locus a 500 kb one would include
    assert classify_location("chr1", 300_000, {}, {}, subtelomere_bp=100_000) != "subtelomeric"
    assert classify_location("chr1", 300_000, {}, {}, subtelomere_bp=500_000) == "subtelomeric"

    # zero reverts to the telomere gap, which no real locus reaches
    assert classify_location("chr1", 10757, {}, {}, subtelomere_bp=0) != "subtelomeric"
    assert classify_location("chr1", 5000, {}, {}, subtelomere_bp=0) == "subtelomeric"


def test_subtelomeric_outranks_a_gene_overlap():
    genes = {"chr1": {"NEAREND": {"start": 20_000, "end": 60_000}}}
    exons = {"chr1": {"NEAREND": {"exon1": {"start": 21_000, "end": 21_200}}}}
    # inside an exon, but within the subtelomeric window, which has priority
    assert classify_location("chr1", 21_100, genes, exons) == "subtelomeric"


def test_chromosome_lengths_agree_with_the_telomere_gaps():
    """The q gap of every chromosome must end exactly at its length."""
    for chrom, length in GRCH38_CHROM_LENGTHS.items():
        q = HG38_TELOMERES[chrom]["q"]
        assert q is not None, chrom
        assert q[1] == length, (chrom, q[1], length)
        assert q[0] == length - 10_000, (chrom, q[0])


def test_an_unknown_contig_is_not_subtelomeric():
    assert classify_location("chrUn_random", 5000, {}, {}) == "intergenic_other"
