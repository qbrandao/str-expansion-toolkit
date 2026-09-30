"""
Checks the orchestration logic (idempotence "skip if already done",
output file names) without running the real bioinformatics tools:
str_toolkit.utils.run_in_env is mocked, and the presence of output
files is simulated.
"""

import shutil
from pathlib import Path

import pytest

from str_toolkit import detect
from str_toolkit.config import Config


def _make_cfg():
    from str_toolkit.config import VamosConfig, TrgtConfig, TandemGenotypesConfig, LongTRConfig
    return Config(
        reference="ref.fa",
        vamos=VamosConfig(
            env_clair3="clair3", env_whatshap="whatshap-env", env_vamos="vamos",
            model_prefix="model", catalog="catalog.tsv",
        ),
        trgt=TrgtConfig(env="trgt", mmi="ref.mmi", repeats_bed="repeats.bed"),
        tandem_genotypes=TandemGenotypesConfig(
            env_last="last_env", env_tandem="tandem-env",
            last_ref_db="lastdb", repeats_bed="repeats.trf.bed",
        ),
        longtr=LongTRConfig(env="longtr", mmi="ref.mmi", regions_bed="longtr.regions.bed"),
    )


def test_run_trgt_skips_if_vcf_already_present(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(detect, "run_in_env", lambda *a, **k: calls.append((a, k)))

    sid = "p01"
    (tmp_path / f"{sid}.sorted.bam").touch()
    (tmp_path / f"{sid}.sorted.bam.bai").touch()
    (tmp_path / f"{sid}.trgt.vcf.gz").touch()

    sample = detect.Sample(sample_id=sid, fastq_path="reads.fastq.gz")
    out = detect.run_trgt(sample, _make_cfg(), tmp_path, threads=4)

    assert out == tmp_path / f"{sid}.trgt.vcf.gz"
    assert calls == []  # nothing re-run, everything already existed


def test_run_trgt_requires_fastq_when_no_bam_cached(tmp_path, monkeypatch):
    monkeypatch.setattr(detect, "run_in_env", lambda *a, **k: None)
    sample = detect.Sample(sample_id="p01")  # neither bam nor fastq
    with pytest.raises(SystemExit):
        detect.run_trgt(sample, _make_cfg(), tmp_path, threads=4)


def test_run_tandem_genotypes_skips_if_tsv_present(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(detect, "run_in_env", lambda *a, **k: calls.append((a, k)))

    sid = "p01"
    (tmp_path / f"{sid}.tandem_genotypes.tsv").touch()

    sample = detect.Sample(sample_id=sid, fastq_path="reads.fastq.gz")
    out = detect.run_tandem_genotypes(sample, _make_cfg(), tmp_path, threads=4)

    assert out == tmp_path / f"{sid}.tandem_genotypes.tsv"
    assert calls == []


def test_run_vamos_requires_bam_or_fastq(tmp_path, monkeypatch):
    monkeypatch.setattr(detect, "run_in_env", lambda *a, **k: None)
    sample = detect.Sample(sample_id="p01")  # neither bam nor fastq
    with pytest.raises(SystemExit):
        detect.run_vamos(sample, _make_cfg(), tmp_path, threads=4)


# ---------------------------------------------------------------------
# Working from an input BAM alone (no FASTQ), which is the common case for
# a BAM received from a sequencing provider.
# ---------------------------------------------------------------------

def _make_indexed_bam(tmp_path: Path, name: str = "input.bam") -> str:
    bam = tmp_path / name
    bam.touch()
    Path(f"{bam}.bai").touch()
    return str(bam)


def test_longtr_reuses_input_bam_instead_of_realigning(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(detect, "run_in_env", lambda env, cmd, **k: calls.append((env, cmd, k)))

    bam = _make_indexed_bam(tmp_path)
    sample = detect.Sample(sample_id="p01", bam_path=bam)  # no fastq at all
    out = detect.run_longtr(sample, _make_cfg(), tmp_path, threads=4)

    assert out == tmp_path / "p01.longtr.vcf.gz"
    joined = " ".join(" ".join(map(str, c[1])) + str(c[2]) for c in calls)
    assert "minimap2" not in joined  # the input BAM was used as is
    assert "LongTR" in joined
    # LongTR was pointed at the input BAM, not at a freshly aligned one
    longtr_call = next(c for c in calls if c[1] and c[1][0] == "LongTR")
    assert bam in longtr_call[1]


def test_realign_forces_minimap2_even_with_an_input_bam(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(detect, "run_in_env", lambda env, cmd, **k: calls.append((env, cmd, k)))

    bam = _make_indexed_bam(tmp_path)
    sample = detect.Sample(sample_id="p01", bam_path=bam, fastq_path="reads.fastq.gz")
    detect.run_longtr(sample, _make_cfg(), tmp_path, threads=4, realign=True)

    pipelines = [c[2].get("shell_pipeline", "") for c in calls]
    assert any("minimap2" in p for p in pipelines)


def test_unindexed_input_bam_is_symlinked_and_indexed(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(detect, "run_in_env", lambda env, cmd, **k: calls.append((env, cmd, k)))

    src = tmp_path / "readonly"
    src.mkdir()
    bam = src / "input.bam"
    bam.touch()  # no .bai next to it

    outdir = tmp_path / "out"
    outdir.mkdir()
    cfg = _make_cfg()
    cfg.samtools_env = "vamos"  # so samtools goes through the mocked run_in_env

    sample = detect.Sample(sample_id="p01", bam_path=str(bam))
    detect.run_longtr(sample, cfg, outdir, threads=4)

    link = outdir / "p01.input.bam"
    assert link.is_symlink()
    assert link.resolve() == bam.resolve()
    # the symlink is what gets indexed, never the input directory
    index_calls = [c for c in calls if c[1][:2] == ["samtools", "index"]]
    assert index_calls and str(link) in index_calls[0][1]


def test_tandem_genotypes_extracts_a_fastq_from_the_bam(tmp_path, monkeypatch):
    calls = []

    def fake_run_in_env(env, cmd, **k):
        calls.append((env, cmd, k))
        # simulate the extraction writing its temporary output
        pipeline = k.get("shell_pipeline", "")
        if "samtools fastq" in pipeline:
            Path(pipeline.split(">")[-1].strip()).touch()

    monkeypatch.setattr(detect, "run_in_env", fake_run_in_env)

    bam = _make_indexed_bam(tmp_path)
    cfg = _make_cfg()
    cfg.samtools_env = "vamos"

    sample = detect.Sample(sample_id="p01", bam_path=bam)  # no fastq
    detect.run_tandem_genotypes(sample, cfg, tmp_path, threads=4)

    extracted = tmp_path / "p01.reads.fastq.gz"
    assert extracted.exists()
    pipelines = [c[2].get("shell_pipeline", "") for c in calls]
    assert any("samtools fastq" in p and "-T MM,ML" in p and "-F 0x900" in p for p in pipelines)
    # compression is parallel when pigz is installed, single-threaded otherwise
    extraction = next(p for p in pipelines if "samtools fastq" in p)
    assert ("pigz" in extraction) == (shutil.which("pigz") is not None)
    # LAST was then given the extracted FASTQ
    last_train = next(c for c in calls if c[1] and c[1][0] == "last-train")
    assert str(extracted) in last_train[1]


def test_extraction_uses_pigz_when_available(tmp_path, monkeypatch):
    calls = []

    def fake_run_in_env(env, cmd, **k):
        calls.append((env, cmd, k))
        pipeline = k.get("shell_pipeline", "")
        if "samtools fastq" in pipeline:
            Path(pipeline.split(">")[-1].strip()).touch()

    monkeypatch.setattr(detect, "run_in_env", fake_run_in_env)
    monkeypatch.setattr(detect.shutil, "which", lambda name: "/usr/bin/pigz" if name == "pigz" else None)

    bam = _make_indexed_bam(tmp_path)
    cfg = _make_cfg()
    cfg.samtools_env = "vamos"

    detect.run_tandem_genotypes(
        detect.Sample(sample_id="p01", bam_path=bam), cfg, tmp_path, threads=8
    )

    extraction = next(c[2]["shell_pipeline"] for c in calls
                      if "samtools fastq" in c[2].get("shell_pipeline", ""))
    assert "pigz -p 8 -c" in extraction
    assert "gzip" not in extraction


def test_supplied_fastq_wins_over_bam_extraction(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(detect, "run_in_env", lambda env, cmd, **k: calls.append((env, cmd, k)))

    bam = _make_indexed_bam(tmp_path)
    sample = detect.Sample(sample_id="p01", bam_path=bam, fastq_path="reads.fastq.gz")
    detect.run_tandem_genotypes(sample, _make_cfg(), tmp_path, threads=4)

    pipelines = [c[2].get("shell_pipeline", "") for c in calls]
    assert not any("samtools fastq" in p for p in pipelines)
    assert not (tmp_path / "p01.reads.fastq.gz").exists()


def test_run_longtr_skips_if_vcf_already_present(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(detect, "run_in_env", lambda *a, **k: calls.append((a, k)))

    sid = "p01"
    (tmp_path / f"{sid}.sorted.bam").touch()
    (tmp_path / f"{sid}.sorted.bam.bai").touch()
    (tmp_path / f"{sid}.longtr.vcf.gz").touch()

    sample = detect.Sample(sample_id=sid, fastq_path="reads.fastq.gz")
    out = detect.run_longtr(sample, _make_cfg(), tmp_path, threads=4)

    assert out == tmp_path / f"{sid}.longtr.vcf.gz"
    assert calls == []  # nothing re-run, everything already existed


def test_run_longtr_requires_fastq_when_no_bam_cached(tmp_path, monkeypatch):
    monkeypatch.setattr(detect, "run_in_env", lambda *a, **k: None)
    sample = detect.Sample(sample_id="p01")  # neither bam nor fastq
    with pytest.raises(SystemExit):
        detect.run_longtr(sample, _make_cfg(), tmp_path, threads=4)


def test_run_trgt_and_run_longtr_reuse_same_alignment(tmp_path, monkeypatch):
    """TRGT and LongTR must reuse the same .sorted.bam if it already exists
    (no duplicate alignment when both run in the same job)."""
    calls = []
    monkeypatch.setattr(detect, "run_in_env", lambda *a, **k: calls.append((a, k)))

    sid = "p01"
    sample = detect.Sample(sample_id=sid, fastq_path="reads.fastq.gz")
    cfg = _make_cfg()

    detect.run_trgt(sample, cfg, tmp_path, threads=4)
    align_calls_after_trgt = sum(1 for a, k in calls if "minimap2" in str(k))
    # run_in_env is mocked (no-op): simulate the .sorted.bam minimap2 would have produced
    (tmp_path / f"{sid}.sorted.bam").touch()

    detect.run_longtr(sample, cfg, tmp_path, threads=4)
    align_calls_after_longtr = sum(1 for a, k in calls if "minimap2" in str(k))

    assert align_calls_after_trgt == 1
    assert align_calls_after_longtr == 1  # no second alignment for LongTR
