"""
Configuration for the detection tools: reference paths, catalogs, and
micromamba environments. Loaded from a YAML file passed via --config.

Example: see config.example.yaml at the repo root.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml


# Every `bin_*` / `bin` field below is the executable to invoke. A bare name
# is resolved inside the matching environment, which is the usual case. An
# absolute path is invoked as is, which is what a tool compiled by hand
# outside any environment needs (vamos and LongTR are commonly built that
# way, e.g. /home/user/vamos/src/vamos). With an absolute path the matching
# `env_*` key can be left empty, and the binary then runs on the current
# PATH rather than through the environment runner.


@dataclass
class VamosConfig:
    env_clair3: str = "clair3"
    env_whatshap: str = "whatshap-env"
    env_vamos: str = "vamos"
    bin_clair3: str = "run_clair3.sh"
    bin_whatshap: str = "whatshap"
    bin_vamos: str = "vamos"
    model_prefix: str = ""
    catalog: str = ""


@dataclass
class TrgtConfig:
    env: str = "trgt"
    bin: str = "trgt"
    mmi: str = ""
    repeats_bed: str = ""


@dataclass
class TandemGenotypesConfig:
    env_last: str = "last_env"
    env_tandem: str = "tandem-env"
    bin_lastal: str = "lastal"
    bin_last_train: str = "last-train"
    bin_last_split: str = "last-split"
    bin_tandem_genotypes: str = "tandem-genotypes"
    last_ref_db: str = ""  # prefix of an index built with `lastdb`
    repeats_bed: str = ""


@dataclass
class LongTRConfig:
    env: str = "longtr"
    bin: str = "LongTR"
    mmi: str = ""  # minimap2 index (can be the same file as trgt.mmi)
    regions_bed: str = ""  # LongTR BED: chrom, start(1-based), end, motif[,motif2], [name]


@dataclass
class Config:
    reference: str = ""
    # Optional environment providing samtools, used for the input inspection
    # and conversion steps that belong to no single tool. Left empty,
    # samtools is taken from the current PATH.
    samtools_env: str = ""
    # How the tool environments are entered: micromamba, mamba, conda, or a
    # full path to one of them. Left empty, the first of the three found on
    # PATH is used.
    env_runner: str = ""
    # Optional environment providing minimap2 AND samtools for the shared
    # alignment step. Left empty, the environment of the tool that asked for
    # the alignment is used, which is the historical behaviour.
    align_env: str = ""
    vamos: VamosConfig = field(default_factory=VamosConfig)
    trgt: TrgtConfig = field(default_factory=TrgtConfig)
    tandem_genotypes: TandemGenotypesConfig = field(default_factory=TandemGenotypesConfig)
    longtr: LongTRConfig = field(default_factory=LongTRConfig)

    @classmethod
    def from_yaml(cls, path: str | Path) -> "Config":
        with open(path) as fh:
            raw = yaml.safe_load(fh) or {}

        return cls(
            reference=raw.get("reference", ""),
            samtools_env=raw.get("samtools_env", ""),
            env_runner=raw.get("env_runner", ""),
            align_env=raw.get("align_env", ""),
            vamos=VamosConfig(**raw.get("vamos", {})),
            trgt=TrgtConfig(**raw.get("trgt", {})),
            tandem_genotypes=TandemGenotypesConfig(**raw.get("tandem_genotypes", {})),
            longtr=LongTRConfig(**raw.get("longtr", {})),
        )
