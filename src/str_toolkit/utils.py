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
    Runs a command inside a given micromamba environment, equivalent to
    `micromamba run -n <env> <cmd...>`.

    - stdout_path: if provided, redirects stdout to this file (equivalent to `> file`).
    - shell_pipeline: if provided (a shell string, e.g. "cmd1 | cmd2 > out"), ignores
      `cmd` and runs this string via `micromamba run -n <env> bash -c "<shell_pipeline>"`.
      Useful for reproducing pipes such as `lastal ... | last-split ...`.
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
    """Runs a command directly on the current PATH (no micromamba run)."""
    logger.info("%s", " ".join(cmd))
    return subprocess.run(cmd, check=check)


def read_samples_list(path: str) -> list[dict]:
    """
    Reads a TSV with a required `sample_id` column, and optional `bam_path` /
    `fastq_path` columns (at least one of the two must be present depending
    on which tools are used). Returns a list of dicts.
    """
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = set(reader.fieldnames or [])
        if "sample_id" not in fieldnames:
            raise SystemExit(f"{path}: missing 'sample_id' column in the TSV")
        if "bam_path" not in fieldnames and "fastq_path" not in fieldnames:
            raise SystemExit(
                f"{path}: at least one of 'bam_path' or 'fastq_path' columns is required"
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
