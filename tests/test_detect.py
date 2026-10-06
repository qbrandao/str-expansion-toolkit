"""
Checks the orchestration logic (idempotence "skip if already done",
output file names) without running the real bioinformatics tools:
str_toolkit.utils.run_in_env is mocked, and the presence of output
files is simulated.
"""

import shutil
import subprocess
from pathlib import Path

import pytest

from str_toolkit import detect
from str_toolkit.config import Config


def _make_cfg():
    from str_toolkit.config import VamosConfig, TrgtConfig, TandemGenotypesConfig, LongTRConfig
    return Config(
        reference="ref.fa",
        # samtools belongs to no single tool, and the indexing and extraction
        # steps route through it, so the fixture names an environment rather
        # than falling through to a real samtools on the PATH
        samtools_env="vamos",
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


def _recording_runner(calls):
    """
    A stand-in for run_in_env that records the call AND creates the file the
    real command would have written.

    detect.py writes every produced file under a temporary name and renames it
    only once the command exits cleanly, so a fake that writes nothing makes
    the rename fail. Creating the target here is what a working tool does.
    """
    def fake(env, cmd=None, **kwargs):
        calls.append((env, cmd if cmd is not None else [], kwargs))
        target = kwargs.get("stdout_path")
        if target:
            Path(target).parent.mkdir(parents=True, exist_ok=True)
            Path(target).touch()
        pipeline = kwargs.get("shell_pipeline") or ""
        out = None
        if ">" in pipeline:
            out = pipeline.rsplit(">", 1)[-1].strip().split()[0]
        elif " -o " in pipeline:            # samtools sort writes with -o
            out = pipeline.rsplit(" -o ", 1)[-1].strip().split()[0]
        if out:
            Path(out).parent.mkdir(parents=True, exist_ok=True)
            Path(out).touch()
    return fake


def test_run_trgt_skips_if_vcf_already_present(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(detect, "run_in_env", _recording_runner(calls))

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
    monkeypatch.setattr(detect, "run_in_env", _recording_runner(calls))

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
    monkeypatch.setattr(detect, "run_in_env", _recording_runner(calls))

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
    monkeypatch.setattr(detect, "run_in_env", _recording_runner(calls))

    bam = _make_indexed_bam(tmp_path)
    sample = detect.Sample(sample_id="p01", bam_path=bam, fastq_path="reads.fastq.gz")
    detect.run_longtr(sample, _make_cfg(), tmp_path, threads=4, realign=True)

    pipelines = [c[2].get("shell_pipeline") or "" for c in calls]
    assert any("minimap2" in p for p in pipelines)


def test_align_env_overrides_the_tool_environment(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(detect, "run_in_env", _recording_runner(calls))

    cfg = _make_cfg()
    cfg.align_env = "bioinfo"   # minimap2 lives here, not in the longtr env
    sample = detect.Sample(sample_id="p01", fastq_path="reads.fastq.gz")
    detect.run_longtr(sample, cfg, tmp_path, threads=4)

    align = next(c for c in calls if "minimap2" in c[2].get("shell_pipeline") or "")
    assert align[0] == "bioinfo"
    # LongTR itself still runs in its own environment
    longtr = next(c for c in calls if c[1] and c[1][0] == "LongTR")
    assert longtr[0] == "longtr"


def test_without_align_env_the_tool_environment_is_used(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(detect, "run_in_env", _recording_runner(calls))

    cfg = _make_cfg()          # align_env left empty
    detect.run_longtr(detect.Sample(sample_id="p01", fastq_path="r.fastq.gz"),
                      cfg, tmp_path, threads=4)

    align = next(c for c in calls if "minimap2" in c[2].get("shell_pipeline") or "")
    assert align[0] == "longtr"


def test_empty_align_env_and_no_tool_env_falls_back_to_path(tmp_path, monkeypatch):
    """A blank env must run on the PATH rather than through the runner."""
    bash_calls = []
    monkeypatch.setattr(detect, "run_in_env",
                        lambda *a, **k: pytest.fail("should not enter an environment"))

    class Done:
        returncode = 0

    monkeypatch.setattr(detect.subprocess, "run",
                        lambda cmd, **k: (bash_calls.append(cmd), Done())[1])
    detect._run_maybe_env("", [], shell_pipeline="minimap2 ... | samtools sort ...")

    assert bash_calls[0][0] == "bash"
    assert "minimap2" in bash_calls[0][2]


def test_unindexed_input_bam_is_symlinked_and_indexed(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(detect, "run_in_env", _recording_runner(calls))

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

    monkeypatch.setattr(detect, "run_in_env", _recording_runner(calls))

    bam = _make_indexed_bam(tmp_path)
    cfg = _make_cfg()
    cfg.samtools_env = "vamos"

    sample = detect.Sample(sample_id="p01", bam_path=bam)  # no fastq
    detect.run_tandem_genotypes(sample, cfg, tmp_path, threads=4)

    extracted = tmp_path / "p01.reads.fastq.gz"
    assert extracted.exists()
    pipelines = [c[2].get("shell_pipeline") or "" for c in calls]
    assert any("samtools fastq" in p and "-T MM,ML" in p and "-F 0x900" in p for p in pipelines)
    # compression is parallel when pigz is installed, single-threaded otherwise
    extraction = next(p for p in pipelines if "samtools fastq" in p)
    assert ("pigz" in extraction) == (shutil.which("pigz") is not None)
    # LAST was then given the extracted FASTQ
    last_train = next(c for c in calls if c[1] and c[1][0] == "last-train")
    assert str(extracted) in last_train[1]


def test_extraction_uses_pigz_when_available(tmp_path, monkeypatch):
    calls = []

    monkeypatch.setattr(detect, "run_in_env", _recording_runner(calls))
    monkeypatch.setattr(detect.shutil, "which", lambda name: "/usr/bin/pigz" if name == "pigz" else None)

    bam = _make_indexed_bam(tmp_path)
    cfg = _make_cfg()
    cfg.samtools_env = "vamos"

    detect.run_tandem_genotypes(
        detect.Sample(sample_id="p01", bam_path=bam), cfg, tmp_path, threads=8
    )

    extraction = next(c[2]["shell_pipeline"] for c in calls
                      if "samtools fastq" in c[2].get("shell_pipeline") or "")
    assert "pigz -p 8 -c" in extraction
    assert "gzip" not in extraction


def test_supplied_fastq_wins_over_bam_extraction(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(detect, "run_in_env", _recording_runner(calls))

    bam = _make_indexed_bam(tmp_path)
    sample = detect.Sample(sample_id="p01", bam_path=bam, fastq_path="reads.fastq.gz")
    detect.run_tandem_genotypes(sample, _make_cfg(), tmp_path, threads=4)

    pipelines = [c[2].get("shell_pipeline") or "" for c in calls]
    assert not any("samtools fastq" in p for p in pipelines)
    assert not (tmp_path / "p01.reads.fastq.gz").exists()


def test_run_longtr_skips_if_vcf_already_present(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(detect, "run_in_env", _recording_runner(calls))

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
    monkeypatch.setattr(detect, "run_in_env", _recording_runner(calls))

    sid = "p01"
    sample = detect.Sample(sample_id=sid, fastq_path="reads.fastq.gz")
    cfg = _make_cfg()

    detect.run_trgt(sample, cfg, tmp_path, threads=4)
    align_calls_after_trgt = sum(1 for _env, _cmd, k in calls if "minimap2" in str(k))
    # the recording runner already created the .sorted.bam, which is what the
    # reuse check looks for
    assert (tmp_path / f"{sid}.sorted.bam").exists()

    detect.run_longtr(sample, cfg, tmp_path, threads=4)
    align_calls_after_longtr = sum(1 for _env, _cmd, k in calls if "minimap2" in str(k))

    assert align_calls_after_trgt == 1
    assert align_calls_after_longtr == 1  # no second alignment for LongTR


# ---------------------------------------------------------------------
# Binaries compiled by hand, outside any environment. vamos and LongTR are
# commonly built that way, e.g. /home/user/vamos/src/vamos.
# ---------------------------------------------------------------------

def test_absolute_vamos_binary_with_no_environment(tmp_path, monkeypatch):
    """An absolute bin path and an empty env must run straight on the PATH."""
    env_calls, bash_calls = [], []

    class Done:
        returncode = 0

    monkeypatch.setattr(detect, "run_in_env",
                        lambda env, cmd, **k: env_calls.append((env, cmd, k)))
    monkeypatch.setattr(detect.subprocess, "run",
                        lambda cmd, **k: (bash_calls.append(cmd), Done())[1])
    monkeypatch.setattr(detect, "run_cmd", lambda cmd, **k: bash_calls.append(cmd))

    cfg = _make_cfg()
    cfg.vamos.env_vamos = ""                                  # no environment
    cfg.vamos.bin_vamos = "/home/user/vamos/src/vamos"        # hand-built

    sid = "p01"
    bam = _make_indexed_bam(tmp_path)
    # everything before the vamos step already exists, so only it runs
    (tmp_path / "phased_merge_output.vcf.gz").touch()
    (tmp_path / f"{sid}_haplotagged.bam").touch()
    (tmp_path / f"{sid}_h1.bam").touch()
    (tmp_path / f"{sid}_h2.bam").touch()

    detect.run_vamos(detect.Sample(sample_id=sid, bam_path=bam), cfg, tmp_path, threads=4)

    # the vamos calls went to the PATH, not through the environment runner
    assert not any(c[1] and str(c[1][0]).endswith("vamos") for c in env_calls)
    vamos_calls = [c for c in bash_calls if c and c[0] == "/home/user/vamos/src/vamos"]
    assert len(vamos_calls) == 2          # one per haplotype
    assert vamos_calls[0][1] == "--contig"


def test_absolute_binary_inside_an_environment_still_uses_the_runner(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(detect, "run_in_env", _recording_runner(calls))

    cfg = _make_cfg()
    cfg.longtr.bin = "/opt/LongTR/LongTR"     # absolute, but env kept

    detect.run_longtr(detect.Sample(sample_id="p01", bam_path=_make_indexed_bam(tmp_path)),
                      cfg, tmp_path, threads=4)

    longtr = next(c for c in calls if c[1] and c[1][0] == "/opt/LongTR/LongTR")
    assert longtr[0] == "longtr"


def test_last_binary_overrides_reach_the_shell_pipeline(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(detect, "run_in_env", _recording_runner(calls))

    cfg = _make_cfg()
    cfg.tandem_genotypes.bin_lastal = "/opt/last/bin/lastal"
    cfg.tandem_genotypes.bin_last_split = "/opt/last/bin/last-split"
    cfg.tandem_genotypes.bin_last_train = "/opt/last/bin/last-train"
    cfg.tandem_genotypes.bin_tandem_genotypes = "/opt/tg/tandem-genotypes"

    detect.run_tandem_genotypes(
        detect.Sample(sample_id="p01", fastq_path="reads.fastq.gz"), cfg, tmp_path, threads=4
    )

    train = next(c for c in calls if c[1] and c[1][0] == "/opt/last/bin/last-train")
    assert "-Q0" in train[1]
    pipeline = next(c[2]["shell_pipeline"] for c in calls
                    if "lastal" in (c[2].get("shell_pipeline") or ""))
    assert "/opt/last/bin/lastal" in pipeline
    assert "/opt/last/bin/last-split" in pipeline
    tg = next(c for c in calls if c[1] and c[1][0] == "/opt/tg/tandem-genotypes")
    assert tg[0] == "tandem-env"


def test_config_reads_binary_overrides(tmp_path):
    from str_toolkit.config import Config

    cfg_file = tmp_path / "config.yaml"
    cfg_file.write_text(
        "reference: ref.fa\n"
        "vamos:\n"
        "  env_vamos: ''\n"
        "  bin_vamos: /home/qbrandao/vamos/src/vamos\n"
        "tandem_genotypes:\n"
        "  env_last: base\n"
        "  env_tandem: base\n"
    )
    cfg = Config.from_yaml(cfg_file)
    assert cfg.vamos.bin_vamos == "/home/qbrandao/vamos/src/vamos"
    assert cfg.vamos.env_vamos == ""
    assert cfg.tandem_genotypes.env_last == "base"
    # untouched keys keep their defaults
    assert cfg.vamos.bin_clair3 == "run_clair3.sh"
    assert cfg.longtr.bin == "LongTR"


# ---------------------------------------------------------------------
# A killed stage in a pipeline must fail the step, and must not leave a
# truncated output behind that the idempotence check would accept.
#
# This is the failure that a real run hit: lastal was killed by the kernel
# for exhausting memory, last-split exited 0 on the short input, the MAF was
# silently truncated, and tandem-genotypes produced a TSV from partial
# alignments that no later run would ever recompute.
# ---------------------------------------------------------------------

def test_wrap_pipeline_is_idempotent():
    from str_toolkit.utils import wrap_pipeline

    assert wrap_pipeline("a | b").startswith("set -o pipefail;")
    already = "set -o pipefail; a | b"
    assert wrap_pipeline(already) == already


def test_a_killed_lastal_leaves_no_maf_and_no_tsv(tmp_path, monkeypatch):
    """The step must fail, and neither the MAF nor the TSV may appear."""
    def dying_runner(env, cmd=None, **kwargs):
        pipeline = kwargs.get("shell_pipeline") or ""
        if "lastal" in pipeline:
            # what pipefail now surfaces instead of a silent success
            raise subprocess.CalledProcessError(137, "lastal")
        target = kwargs.get("stdout_path")
        if target:
            Path(target).touch()

    monkeypatch.setattr(detect, "run_in_env", dying_runner)

    sid = "p01"
    with pytest.raises(subprocess.CalledProcessError):
        detect.run_tandem_genotypes(
            detect.Sample(sample_id=sid, fastq_path="reads.fastq.gz"),
            _make_cfg(), tmp_path, threads=4,
        )

    assert not (tmp_path / f"{sid}_alignments.maf").exists()
    assert not (tmp_path / f"{sid}.tandem_genotypes.tsv").exists()
    # the parameters file finished before the failure and is kept, so a rerun
    # does not repeat last-train
    assert (tmp_path / f"{sid}_reads.par").exists()


def test_a_completed_maf_is_not_recomputed(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(detect, "run_in_env", _recording_runner(calls))

    sid = "p01"
    (tmp_path / f"{sid}_reads.par").touch()
    (tmp_path / f"{sid}_alignments.maf").touch()

    detect.run_tandem_genotypes(
        detect.Sample(sample_id=sid, fastq_path="reads.fastq.gz"),
        _make_cfg(), tmp_path, threads=4,
    )

    pipelines = [c[2].get("shell_pipeline") or "" for c in calls]
    assert not any("lastal" in p for p in pipelines)
    assert not any(c[1] and "last-train" in str(c[1][0]) for c in calls)
    # only the genotyping step ran
    assert any(c[1] and c[1][0] == "tandem-genotypes" for c in calls)


def test_last_memory_knobs_reach_the_command(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(detect, "run_in_env", _recording_runner(calls))

    cfg = _make_cfg()
    cfg.tandem_genotypes.last_threads = 6        # fewer than --threads
    cfg.tandem_genotypes.lastal_batch_size = "1G"
    cfg.tandem_genotypes.lastal_extra_args = "-C2"

    detect.run_tandem_genotypes(
        detect.Sample(sample_id="p01", fastq_path="reads.fastq.gz"),
        cfg, tmp_path, threads=16,
    )

    pipeline = next(c[2]["shell_pipeline"] for c in calls
                    if "lastal" in (c[2].get("shell_pipeline") or ""))
    assert "-P6" in pipeline and "-P16" not in pipeline
    assert "-i 1G" in pipeline
    assert "-C2" in pipeline

    train = next(c for c in calls if c[1] and c[1][0] == "last-train")
    assert "6" in train[1]        # last-train uses the same reduced count


def test_last_threads_defaults_to_the_run_threads(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(detect, "run_in_env", _recording_runner(calls))

    detect.run_tandem_genotypes(
        detect.Sample(sample_id="p01", fastq_path="reads.fastq.gz"),
        _make_cfg(), tmp_path, threads=16,
    )

    pipeline = next(c[2]["shell_pipeline"] for c in calls
                    if "lastal" in (c[2].get("shell_pipeline") or ""))
    assert "-P16" in pipeline
    assert "-i " not in pipeline          # no batch bound unless configured


def test_alignment_bam_is_atomic(tmp_path, monkeypatch):
    """A killed minimap2 must not leave a .sorted.bam that looks finished."""
    def dying_runner(env, cmd=None, **kwargs):
        if "minimap2" in (kwargs.get("shell_pipeline") or ""):
            raise subprocess.CalledProcessError(137, "minimap2")

    monkeypatch.setattr(detect, "run_in_env", dying_runner)

    with pytest.raises(subprocess.CalledProcessError):
        detect.run_longtr(detect.Sample(sample_id="p01", fastq_path="r.fastq.gz"),
                          _make_cfg(), tmp_path, threads=4)

    assert not (tmp_path / "p01.sorted.bam").exists()


def test_vamos_short_circuits_on_finished_haplotype_vcfs(tmp_path, monkeypatch):
    """
    A finished VAMOS run must not rebuild its intermediates.

    The haplotagged BAM and the two haplotype BAMs together weigh about twice
    the input BAM, so they get deleted after a run. The per-step skips are
    sequential, so without a short circuit their absence would restart
    whatshap even though the VCFs they lead to exist.
    """
    calls = []
    monkeypatch.setattr(detect, "run_in_env", _recording_runner(calls))

    sid = "p01"
    (tmp_path / f"{sid}_assembly.hap1.vcf").touch()
    (tmp_path / f"{sid}_assembly.hap2.vcf").touch()
    # every intermediate deleted, as a cleanup would leave it
    assert not (tmp_path / f"{sid}_haplotagged.bam").exists()

    out = detect.run_vamos(detect.Sample(sample_id=sid, bam_path=_make_indexed_bam(tmp_path)),
                           _make_cfg(), tmp_path, threads=4)

    assert out == {"hap1": tmp_path / f"{sid}_assembly.hap1.vcf",
                   "hap2": tmp_path / f"{sid}_assembly.hap2.vcf"}
    assert calls == []


def test_vamos_still_runs_when_only_one_haplotype_vcf_exists(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(detect, "run_in_env", _recording_runner(calls))

    sid = "p01"
    (tmp_path / f"{sid}_assembly.hap1.vcf").touch()   # hap2 missing
    (tmp_path / "phased_merge_output.vcf.gz").touch()
    (tmp_path / f"{sid}_haplotagged.bam").touch()
    (tmp_path / f"{sid}_h1.bam").touch()
    (tmp_path / f"{sid}_h2.bam").touch()

    detect.run_vamos(detect.Sample(sample_id=sid, bam_path=_make_indexed_bam(tmp_path)),
                     _make_cfg(), tmp_path, threads=4)

    # only the missing haplotype was genotyped
    vamos_calls = [c for c in calls if c[1] and c[1][0] == "vamos"]
    assert len(vamos_calls) == 1
    assert str(tmp_path / f"{sid}_assembly.hap2.vcf") in vamos_calls[0][1]
