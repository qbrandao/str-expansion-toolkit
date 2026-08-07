from str_toolkit.cli import build_parser


def test_detect_requires_sample_or_list():
    parser = build_parser()
    try:
        parser.parse_args(["detect", "--config", "config.yaml", "-o", "out/"])
        assert False, "should fail without --sample or --samples-list"
    except SystemExit:
        pass


def test_detect_single_sample_parses():
    parser = build_parser()
    args = parser.parse_args(
        [
            "detect",
            "--sample", "p01",
            "--bam", "p01.bam",
            "--fastq", "p01.fastq.gz",
            "--config", "config.yaml",
            "-o", "out/",
        ]
    )
    assert args.sample == "p01"
    assert args.tools == ["vamos", "tandem-genotypes", "longtr"]


def test_detect_trgt_available_as_explicit_opt_in():
    parser = build_parser()
    args = parser.parse_args(
        [
            "detect",
            "--sample", "p01",
            "--bam", "p01.bam",
            "--fastq", "p01.fastq.gz",
            "--config", "config.yaml",
            "-o", "out/",
            "--tools", "vamos", "trgt", "tandem-genotypes", "longtr",
        ]
    )
    assert "trgt" in args.tools


def test_compare_parses():
    parser = build_parser()
    args = parser.parse_args(
        [
            "compare",
            "--patients-dir", "results/patients",
            "--controls-json", "controls.json",
            "--genes-bed", "genes.bed.gz",
            "--exons-bed", "exons.bed.gz",
            "-o", "report.tsv",
        ]
    )
    assert args.patients_dir == "results/patients"
    assert args.format == "tsv"
    assert args.threshold == 0
