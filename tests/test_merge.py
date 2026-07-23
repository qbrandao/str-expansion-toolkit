from str_toolkit.merge import STRCall, build_locus_record, canonical_motif, cluster_calls


def test_canonical_motif_rotation_equivalence():
    assert canonical_motif("AAAG") == canonical_motif("AAGA") == canonical_motif("AGAA") == canonical_motif("GAAA")


def test_canonical_motif_case_insensitive():
    assert canonical_motif("aaag") == canonical_motif("AAAG")


def test_cluster_merges_shifted_coordinates_same_motif():
    # Même locus biologique, ancré différemment par VAMOS/TRGT/tandem-genotypes
    calls = [
        STRCall(chrom="chr1", start=1000, end=1040, motif="AAAG", size=42, source="vamos_hap1"),
        STRCall(chrom="chr1", start=1006, end=1046, motif="AAGA", size=38, source="trgt_allele1"),
        STRCall(chrom="chr1", start=1012, end=1050, motif="GAAA", size=10, source="tandem_genotypes_allele1"),
    ]
    clusters = cluster_calls(calls, window=25)
    assert len(clusters) == 1
    assert len(clusters[0]) == 3


def test_cluster_does_not_merge_different_motifs_at_same_position():
    calls = [
        STRCall(chrom="chr1", start=1000, end=1040, motif="AAAG", size=42, source="vamos_hap1"),
        STRCall(chrom="chr1", start=1005, end=1045, motif="AC", size=38, source="trgt_allele1"),
    ]
    clusters = cluster_calls(calls, window=25)
    assert len(clusters) == 2


def test_cluster_does_not_merge_distant_loci():
    calls = [
        STRCall(chrom="chr1", start=1000, end=1040, motif="AAAG", size=42, source="vamos_hap1"),
        STRCall(chrom="chr1", start=50000, end=50040, motif="AAAG", size=38, source="trgt_allele1"),
    ]
    clusters = cluster_calls(calls, window=25)
    assert len(clusters) == 2


def test_cluster_respects_chrom_boundaries():
    calls = [
        STRCall(chrom="chr1", start=1000, end=1040, motif="AAAG", size=42, source="vamos_hap1"),
        STRCall(chrom="chr2", start=1000, end=1040, motif="AAAG", size=38, source="trgt_allele1"),
    ]
    clusters = cluster_calls(calls, window=25)
    assert len(clusters) == 2


def test_build_locus_record_keeps_sizes_separate_by_source():
    cluster = [
        STRCall(chrom="chr1", start=1000, end=1040, motif="AAAG", size=42, source="vamos_hap1"),
        STRCall(chrom="chr1", start=1006, end=1046, motif="AAGA", size=38, source="trgt_allele1"),
    ]
    record = build_locus_record(cluster)
    assert record["chrom"] == "chr1"
    assert record["pos"] == 1000
    assert record["end"] == 1046
    assert record["sizes_by_source"] == {"vamos_hap1": 42, "trgt_allele1": 38}
    assert set(record["sources"]) == {"vamos_hap1", "trgt_allele1"}
