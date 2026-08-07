"""
Configuration for the detection tools: reference paths, catalogs, and
micromamba environments. Loaded from a YAML file passed via --config.

Example: see config.example.yaml at the repo root.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml


@dataclass
class VamosConfig:
    env_clair3: str = "clair3"
    env_whatshap: str = "whatshap-env"
    env_vamos: str = "vamos"
    model_prefix: str = ""
    catalog: str = ""


@dataclass
class TrgtConfig:
    env: str = "trgt"
    mmi: str = ""
    repeats_bed: str = ""


@dataclass
class TandemGenotypesConfig:
    env_last: str = "last_env"
    env_tandem: str = "tandem-env"
    last_ref_db: str = ""  # prefix of an index built with `lastdb`
    repeats_bed: str = ""


@dataclass
class LongTRConfig:
    env: str = "longtr"
    mmi: str = ""  # minimap2 index (can be the same file as trgt.mmi)
    regions_bed: str = ""  # LongTR BED: chrom, start(1-based), end, motif[,motif2], [name]


@dataclass
class Config:
    reference: str = ""
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
            vamos=VamosConfig(**raw.get("vamos", {})),
            trgt=TrgtConfig(**raw.get("trgt", {})),
            tandem_genotypes=TandemGenotypesConfig(**raw.get("tandem_genotypes", {})),
            longtr=LongTRConfig(**raw.get("longtr", {})),
        )
