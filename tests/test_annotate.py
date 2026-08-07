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
