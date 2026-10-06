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


def test_meiotic_instability_excludes_sex_chromosomes_by_default():
    parser = build_parser()
    args = parser.parse_args([
        "meiotic-instability",
        "--duos", "duos.tsv", "--data-dir", "results/",
        "--genes-bed", "genes.bed.gz", "--exons-bed", "exons.bed.gz",
        "-o", "out.tsv",
    ])
    assert args.include_sex_chromosomes is False


def test_meiotic_instability_include_sex_chromosomes_flag():
    parser = build_parser()
    args = parser.parse_args([
        "meiotic-instability",
        "--duos", "duos.tsv", "--data-dir", "results/",
        "--genes-bed", "genes.bed.gz", "--exons-bed", "exons.bed.gz",
        "-o", "out.tsv", "--include-sex-chromosomes",
    ])
    assert args.include_sex_chromosomes is True


# ---------------------------------------------------------------------
# --version, which exists to tell a stale install from a current one.
# ---------------------------------------------------------------------

def test_version_reports_the_module_path_and_features(capsys):
    from str_toolkit.cli import main

    assert main(["--version"]) == 0
    out = capsys.readouterr().out
    assert "str-expansion-toolkit" in out
    assert "module" in out and "str_toolkit" in out
    assert "interpreter" in out
    # the features an old install would be missing
    assert "repertoire --subtelomere-bp" in out
    assert "detect --realign" in out


def test_version_needs_no_subcommand(capsys):
    from str_toolkit.cli import main

    # --version alone must not fail on a missing subcommand
    assert main(["--version"]) == 0


def test_no_subcommand_prints_help_and_fails(capsys):
    from str_toolkit.cli import main

    assert main([]) == 2
    assert "usage: str-toolkit" in capsys.readouterr().out


def test_repertoire_accepts_subtelomere_bp():
    """The option whose absence signalled a stale install."""
    from str_toolkit.cli import build_parser

    args = build_parser().parse_args([
        "repertoire", "--controls-dir", "d", "--genes-bed", "g.gz",
        "--exons-bed", "e.gz", "-o", "out.tsv", "--subtelomere-bp", "250000",
    ])
    assert args.subtelomere_bp == 250_000


def test_instability_subcommands_accept_subtelomere_bp():
    from str_toolkit.cli import build_parser

    parser = build_parser()
    meiotic = parser.parse_args([
        "meiotic-instability", "--duos", "d.tsv", "--data-dir", "x",
        "--genes-bed", "g.gz", "--exons-bed", "e.gz", "-o", "m.tsv",
        "--subtelomere-bp", "100000",
    ])
    assert meiotic.subtelomere_bp == 100_000

    somatic = parser.parse_args([
        "somatic-instability", "--samples-list", "s.tsv", "--detect-dir", "x",
        "--genes-bed", "g.gz", "--exons-bed", "e.gz", "-o", "s.tsv",
        "--subtelomere-bp", "5000000",
    ])
    assert somatic.subtelomere_bp == 5_000_000
