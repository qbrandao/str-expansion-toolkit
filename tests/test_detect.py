"""
Vérifie la logique d'orchestration (idempotence "skip si déjà fait",
noms de fichiers produits) sans exécuter les vrais outils bioinfo :
on mocke str_toolkit.utils.run_in_env et on simule la présence de
fichiers de sortie.
"""

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
    assert calls == []  # rien relancé, tout existait déjà


def test_run_trgt_requires_fastq_when_no_bam_cached(tmp_path, monkeypatch):
    monkeypatch.setattr(detect, "run_in_env", lambda *a, **k: None)
    sample = detect.Sample(sample_id="p01")  # ni bam ni fastq
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


def test_run_vamos_requires_bam(tmp_path, monkeypatch):
    monkeypatch.setattr(detect, "run_in_env", lambda *a, **k: None)
    sample = detect.Sample(sample_id="p01", fastq_path="reads.fastq.gz")  # pas de bam
    with pytest.raises(SystemExit):
        detect.run_vamos(sample, _make_cfg(), tmp_path, threads=4)


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
    assert calls == []  # rien relancé, tout existait déjà


def test_run_longtr_requires_fastq_when_no_bam_cached(tmp_path, monkeypatch):
    monkeypatch.setattr(detect, "run_in_env", lambda *a, **k: None)
    sample = detect.Sample(sample_id="p01")  # ni bam ni fastq
    with pytest.raises(SystemExit):
        detect.run_longtr(sample, _make_cfg(), tmp_path, threads=4)


def test_run_trgt_and_run_longtr_reuse_same_alignment(tmp_path, monkeypatch):
    """TRGT et LongTR doivent réutiliser le même .sorted.bam s'il existe déjà
    (pas de ré-alignement en double quand les deux tournent dans le même run)."""
    calls = []
    monkeypatch.setattr(detect, "run_in_env", lambda *a, **k: calls.append((a, k)))

    sid = "p01"
    sample = detect.Sample(sample_id=sid, fastq_path="reads.fastq.gz")
    cfg = _make_cfg()

    detect.run_trgt(sample, cfg, tmp_path, threads=4)
    align_calls_after_trgt = sum(1 for a, k in calls if "minimap2" in str(k))
    # run_in_env est mocké (no-op) : simule le .sorted.bam que minimap2 aurait produit
    (tmp_path / f"{sid}.sorted.bam").touch()

    detect.run_longtr(sample, cfg, tmp_path, threads=4)
    align_calls_after_longtr = sum(1 for a, k in calls if "minimap2" in str(k))

    assert align_calls_after_trgt == 1
    assert align_calls_after_longtr == 1  # pas de second alignement pour LongTR
