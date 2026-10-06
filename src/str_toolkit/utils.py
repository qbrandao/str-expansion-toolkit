from __future__ import annotations

import csv
import logging
import shutil
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)


def ensure_outdir(path: Path) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------
# Environment runner
#
# The tool environments are entered through micromamba, mamba or conda,
# depending on what the machine has. All three accept `run -n <env>`, but
# conda buffers the child's output unless told not to, which would hide a
# long tool's progress and has historically masked exit codes, so conda
# gets --no-capture-output.
#
# The runner is chosen once, from the `env_runner` key in config.yaml, and
# falls back to whichever of the three is on PATH when that key is absent.
# ---------------------------------------------------------------------

_SUPPORTED_RUNNERS = ("micromamba", "mamba", "conda")
_env_runner: str | None = None


def set_env_runner(runner: str | None = None) -> str:
    """
    Fixes the runner used by run_in_env. Given a name, it is used as is, so
    an absolute path to a conda that is not on PATH works too. Given
    nothing, the first of micromamba, mamba or conda found on PATH wins.
    Raises SystemExit when none is available.
    """
    global _env_runner

    if runner:
        _env_runner = runner
        if not shutil.which(runner) and not Path(runner).exists():
            logger.warning(
                "env_runner '%s' is neither on PATH nor an existing file. Every tool "
                "invocation will fail until it is installed or the path is corrected.",
                runner,
            )
        return _env_runner

    for candidate in _SUPPORTED_RUNNERS:
        if shutil.which(candidate):
            _env_runner = candidate
            logger.info("environment runner: %s (auto-detected)", candidate)
            return _env_runner

    raise SystemExit(
        "No environment runner found. Install one of "
        f"{', '.join(_SUPPORTED_RUNNERS)}, or set 'env_runner' in config.yaml to "
        "its full path."
    )


def get_env_runner() -> str:
    """Returns the current runner, auto-detecting one on first use."""
    if _env_runner is None:
        return set_env_runner()
    return _env_runner


def _env_prefix(env: str) -> list[str]:
    runner = get_env_runner()
    # basename, so that a full path such as /opt/conda/bin/conda is still
    # recognised as conda and gets --no-capture-output
    if Path(runner).name.startswith("conda"):
        return [runner, "run", "--no-capture-output", "-n", env]
    return [runner, "run", "-n", env]


def run_in_env(
    env: str,
    cmd: list[str],
    *,
    stdout_path: Path | str | None = None,
    shell_pipeline: str | None = None,
    check: bool = True,
) -> subprocess.CompletedProcess:
    """
    Runs a command inside a given conda-style environment, equivalent to
    `<runner> run -n <env> <cmd...>` (see set_env_runner).

    - stdout_path: if provided, redirects stdout to this file (equivalent to `> file`).
    - shell_pipeline: if provided (a shell string, e.g. "cmd1 | cmd2 > out"), ignores
      `cmd` and runs this string via `<runner> run -n <env> bash -c "<shell_pipeline>"`.
      Useful for reproducing pipes such as `lastal ... | last-split ...`.
    """
    prefix = _env_prefix(env)

    if shell_pipeline is not None:
        full_cmd = [*prefix, "bash", "-c", shell_pipeline]
        logger.info("[%s] bash -c: %s", env, shell_pipeline)
        return subprocess.run(full_cmd, check=check)

    full_cmd = [*prefix, *cmd]
    logger.info("[%s] %s", env, " ".join(cmd))

    if stdout_path is not None:
        with open(stdout_path, "w") as out_fh:
            return subprocess.run(full_cmd, stdout=out_fh, check=check)

    return subprocess.run(full_cmd, check=check)


def run_cmd(cmd: list[str], *, check: bool = True) -> subprocess.CompletedProcess:
    """Runs a command directly on the current PATH (no micromamba run)."""
    logger.info("%s", " ".join(cmd))
    return subprocess.run(cmd, check=check)


def read_tsv_dicts(path: str, required_columns: set[str]) -> list[dict]:
    """Generic TSV reader: validates required columns are present, returns rows as dicts."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = set(reader.fieldnames or [])
        missing = required_columns - fieldnames
        if missing:
            raise SystemExit(f"{path}: missing columns: {', '.join(sorted(missing))}")
        return [dict(row) for row in reader]


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
