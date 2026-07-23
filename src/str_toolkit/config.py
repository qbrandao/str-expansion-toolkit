"""
Configuration des 3 outils : chemins de référence, catalogues, environnements
micromamba. Chargée depuis un fichier YAML passé via --config.

Exemple : voir config.example.yaml à la racine du repo.
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
    last_ref_db: str = ""  # préfixe d'un index construit avec `lastdb`
    repeats_bed: str = ""


@dataclass
class Config:
    reference: str = ""
    vamos: VamosConfig = field(default_factory=VamosConfig)
    trgt: TrgtConfig = field(default_factory=TrgtConfig)
    tandem_genotypes: TandemGenotypesConfig = field(default_factory=TandemGenotypesConfig)

    @classmethod
    def from_yaml(cls, path: str | Path) -> "Config":
        with open(path) as fh:
            raw = yaml.safe_load(fh) or {}

        return cls(
            reference=raw.get("reference", ""),
            vamos=VamosConfig(**raw.get("vamos", {})),
            trgt=TrgtConfig(**raw.get("trgt", {})),
            tandem_genotypes=TandemGenotypesConfig(**raw.get("tandem_genotypes", {})),
        )
