from __future__ import annotations

import csv
import logging
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)


def ensure_outdir(path: Path) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def run_in_env(
    env: str,
    cmd: list[str],
    *,
    stdout_path: Path | str | None = None,
    shell_pipeline: str | None = None,
    check: bool = True,
) -> subprocess.CompletedProcess:
    """
    Exécute une commande dans un environnement micromamba donné, équivalent à
    `micromamba run -n <env> <cmd...>`.

    - stdout_path : si fourni, redirige stdout vers ce fichier (équivalent `> file`).
    - shell_pipeline : si fourni (chaîne shell, ex: "cmd1 | cmd2 > out"), ignore `cmd`
      et exécute cette chaîne via `micromamba run -n <env> bash -c "<shell_pipeline>"`.
      Utile pour reproduire des pipes comme `lastal ... | last-split ...`.
    """
    if shell_pipeline is not None:
        full_cmd = ["micromamba", "run", "-n", env, "bash", "-c", shell_pipeline]
        logger.info("[%s] bash -c: %s", env, shell_pipeline)
        return subprocess.run(full_cmd, check=check)

    full_cmd = ["micromamba", "run", "-n", env, *cmd]
    logger.info("[%s] %s", env, " ".join(full_cmd[4:]))

    if stdout_path is not None:
        with open(stdout_path, "w") as out_fh:
            return subprocess.run(full_cmd, stdout=out_fh, check=check)

    return subprocess.run(full_cmd, check=check)


def run_cmd(cmd: list[str], *, check: bool = True) -> subprocess.CompletedProcess:
    """Exécute une commande directement sur le PATH courant (pas de micromamba run)."""
    logger.info("%s", " ".join(cmd))
    return subprocess.run(cmd, check=check)


def read_samples_list(path: str) -> list[dict]:
    """
    Lit un TSV avec la colonne `sample_id` obligatoire, et `bam_path` /
    `fastq_path` optionnelles (au moins l'une des deux doit être présente
    selon les outils utilisés). Retourne une liste de dicts.
    """
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = set(reader.fieldnames or [])
        if "sample_id" not in fieldnames:
            raise SystemExit(f"{path}: colonne 'sample_id' manquante dans le TSV")
        if "bam_path" not in fieldnames and "fastq_path" not in fieldnames:
            raise SystemExit(
                f"{path}: au moins une des colonnes 'bam_path' ou 'fastq_path' est requise"
            )
        for row in reader:
            rows.append(
                {
                    "sample_id": row["sample_id"],
                    "bam_path": row.get("bam_path") or None,
                    "fastq_path": row.get("fastq_path") or None,
                }
            )
    return rows
