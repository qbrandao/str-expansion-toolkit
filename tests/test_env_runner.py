"""
Checks how the tool environments are entered: micromamba, mamba or conda,
chosen from config.yaml or auto-detected. The tools themselves are never
run, only the argv that would be handed to subprocess is inspected.
"""

from __future__ import annotations

import pytest

from str_toolkit import utils


@pytest.fixture(autouse=True)
def _reset_runner():
    """Each test starts with no runner fixed, since it is module state."""
    utils._env_runner = None
    yield
    utils._env_runner = None


def _captured_argv(monkeypatch):
    calls = []

    class Done:
        returncode = 0

    monkeypatch.setattr(utils.subprocess, "run",
                        lambda cmd, **k: (calls.append(cmd), Done())[1])
    return calls


def test_explicit_micromamba(monkeypatch):
    calls = _captured_argv(monkeypatch)
    utils.set_env_runner("micromamba")
    utils.run_in_env("longtr", ["LongTR", "--bams", "x.bam"])

    assert calls[0][:4] == ["micromamba", "run", "-n", "longtr"]
    assert "--no-capture-output" not in calls[0]
    assert calls[0][4:] == ["LongTR", "--bams", "x.bam"]


def test_conda_gets_no_capture_output(monkeypatch):
    calls = _captured_argv(monkeypatch)
    utils.set_env_runner("conda")
    utils.run_in_env("clair3", ["run_clair3.sh", "--platform=ont"])

    assert calls[0][:5] == ["conda", "run", "--no-capture-output", "-n", "clair3"]
    assert calls[0][5:] == ["run_clair3.sh", "--platform=ont"]


def test_full_path_to_conda_is_still_recognised_as_conda(monkeypatch, tmp_path):
    calls = _captured_argv(monkeypatch)
    runner = tmp_path / "anaconda3" / "bin" / "conda"
    runner.parent.mkdir(parents=True)
    runner.touch()

    utils.set_env_runner(str(runner))
    utils.run_in_env("vamos", ["vamos", "--contig"])

    assert calls[0][:3] == [str(runner), "run", "--no-capture-output"]


def test_mamba_is_treated_like_micromamba(monkeypatch):
    calls = _captured_argv(monkeypatch)
    utils.set_env_runner("mamba")
    utils.run_in_env("vamos", ["vamos"])

    assert "--no-capture-output" not in calls[0]
    assert calls[0][:4] == ["mamba", "run", "-n", "vamos"]


def test_autodetection_prefers_micromamba(monkeypatch):
    monkeypatch.setattr(utils.shutil, "which", lambda n: f"/usr/bin/{n}" if n in ("micromamba", "conda") else None)
    assert utils.set_env_runner() == "micromamba"


def test_autodetection_falls_back_to_conda(monkeypatch):
    monkeypatch.setattr(utils.shutil, "which", lambda n: "/opt/conda/bin/conda" if n == "conda" else None)
    assert utils.set_env_runner() == "conda"


def test_no_runner_available_is_a_clear_error(monkeypatch):
    monkeypatch.setattr(utils.shutil, "which", lambda n: None)
    with pytest.raises(SystemExit) as exc:
        utils.set_env_runner()
    assert "env_runner" in str(exc.value)


def test_unknown_runner_warns_but_does_not_abort(monkeypatch, caplog):
    monkeypatch.setattr(utils.shutil, "which", lambda n: None)
    with caplog.at_level("WARNING"):
        assert utils.set_env_runner("nosuchrunner") == "nosuchrunner"
    assert "neither on PATH" in caplog.text


def test_shell_pipeline_keeps_the_prefix(monkeypatch):
    calls = _captured_argv(monkeypatch)
    utils.set_env_runner("conda")
    utils.run_in_env("last_env", [], shell_pipeline="lastal ... | last-split -m1 > out.maf")

    assert calls[0][:5] == ["conda", "run", "--no-capture-output", "-n", "last_env"]
    assert calls[0][5] == "bash"
    assert calls[0][6] == "-c"
    assert "last-split" in calls[0][7]
    # pipefail, so a stage killed mid-pipeline fails the step instead of
    # leaving a truncated output behind a zero exit code
    assert calls[0][7].startswith("set -o pipefail;")


def test_every_pipeline_is_wrapped_once(monkeypatch):
    calls = _captured_argv(monkeypatch)
    utils.set_env_runner("micromamba")
    utils.run_in_env("e", [], shell_pipeline="set -o pipefail; a | b")

    assert calls[0][-1].count("pipefail") == 1


def test_config_env_runner_is_applied(tmp_path):
    from str_toolkit.config import Config

    cfg_file = tmp_path / "config.yaml"
    cfg_file.write_text(
        "reference: ref.fa\n"
        "env_runner: /home/user/anaconda3/bin/conda\n"
        "samtools_env: vamos\n"
    )
    cfg = Config.from_yaml(cfg_file)
    assert cfg.env_runner == "/home/user/anaconda3/bin/conda"
    assert cfg.samtools_env == "vamos"


def test_config_without_env_runner_leaves_it_empty(tmp_path):
    from str_toolkit.config import Config

    cfg_file = tmp_path / "config.yaml"
    cfg_file.write_text("reference: ref.fa\n")
    assert Config.from_yaml(cfg_file).env_runner == ""
