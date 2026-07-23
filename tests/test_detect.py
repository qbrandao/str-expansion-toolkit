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
    from str_toolkit.config import VamosConfig, TrgtConfig, TandemGenotypesConfig
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
