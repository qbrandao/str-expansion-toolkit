"""
Fusion des sorties VAMOS / TRGT / tandem-genotypes en un VCF unique par
échantillon.

PROBLÈME ADRESSÉ : chaque outil ancre ses coordonnées différemment pour ce
qui est biologiquement le même locus STR (VAMOS = coord. d'assemblage,
TRGT = coord. du catalogue bed, tandem-genotypes = coord. TRF). On ne peut
donc pas fusionner par égalité stricte de (chrom, pos, motif) : il faut un
matching tolérant par intervalle + motif canonique.

IMPORTANT sur les unités : les tailles ne sont PAS forcément comparables
entre outils (VAMOS : longueur en unités de motif ; TRGT et tandem-genotypes :
longueur en bp, mais mesurées différemment -- TRGT par génotypage direct,
tandem-genotypes par clustering de longueurs read-level). Ce module ne
fusionne donc PAS les tailles en une valeur unique : chaque locus fusionné
garde une taille par source (SIZES=vamos_hap1:42,trgt_allele1:38,...).
L'interprétation/comparaison aux contrôles se fait en aval, source par
source (voir str_toolkit/compare.py).
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median

import pysam

# ---------------------------------------------------------------------
# Représentation normalisée d'un call STR, indépendante de l'outil source
# ---------------------------------------------------------------------

@dataclass
class STRCall:
    chrom: str
    start: int
    end: int
    motif: str
    size: float
    source: str  # ex: "vamos_hap1", "trgt_allele2", "tandem_genotypes_allele1"
    raw: dict = field(default_factory=dict)


def canonical_motif(motif: str) -> str:
    """Rotation circulaire minimale : AAGA / AGAA / GAAA / AAAG -> AAAG."""
    m = motif.strip().upper()
    if not m:
        return m
    rotations = [m[i:] + m[:i] for i in range(len(m))]
    return min(rotations)


# ---------------------------------------------------------------------
# Parsers par outil — ajuste les noms de champs si tes VCF/TSV diffèrent
# ---------------------------------------------------------------------

def _first(value):
    """INFO peut être une valeur unique ou un tuple selon Number= dans le header."""
    if isinstance(value, (tuple, list)):
        return value[0] if value else None
    return value


def parse_vamos(hap_vcfs: dict[str, Path]) -> list[STRCall]:
    """
    Lit les VCF VAMOS par haplotype (ex: {"hap1": ..., "hap2": ...}).
    Champs utilisés (identiques à STRlist2json.py) : INFO/RU (motif,
    premier élément si liste), INFO/LEN_H1 (taille, en unités de motif).
    """
    calls = []
    for hap, vcf_path in hap_vcfs.items():
        vcf_path = Path(vcf_path)
        if not vcf_path.exists():
            continue
        with pysam.VariantFile(str(vcf_path)) as vf:
            for record in vf:
                ru = record.info.get("RU")
                length = record.info.get("LEN_H1")
                if ru is None or length is None:
                    continue
                motif = str(_first(ru)).split(",")[0]
                end = record.info.get("END", record.stop)
                calls.append(
                    STRCall(
                        chrom=record.chrom,
                        start=record.pos,
                        end=int(end) if end else record.pos,
                        motif=motif,
                        size=float(_first(length)),
                        source=f"vamos_{hap}",
                        raw={"info": dict(record.info)},
                    )
                )
    return calls


def parse_trgt(vcf_path: Path) -> list[STRCall]:
    """
    Lit le VCF produit par `trgt genotype` (.vcf.gz).

    ATTENTION : les noms de champs ci-dessous (INFO/MOTIFS ou INFO/RU,
    FORMAT/AL) sont ceux généralement utilisés par TRGT, mais peuvent varier
    selon la version. Vérifie avec `zcat sample.trgt.vcf.gz | grep '^##'`
    et ajuste si besoin.
    """
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        return []
    calls = []
    with pysam.VariantFile(str(vcf_path)) as vf:
        for record in vf:
            motif_field = record.info.get("MOTIFS") or record.info.get("RU")
            if motif_field is None:
                continue
            motif = str(_first(motif_field)).split(",")[0]
            end = record.info.get("END", record.stop)

            for sample_name, sample_data in record.samples.items():
                al = sample_data.get("AL")  # allele lengths en bp, ex: (10, 14)
                if al is None:
                    continue
                for i, allele_len in enumerate(al):
                    if allele_len is None:
                        continue
                    calls.append(
                        STRCall(
                            chrom=record.chrom,
                            start=record.pos,
                            end=int(end) if end else record.pos,
                            motif=motif,
                            size=float(allele_len),
                            source=f"trgt_allele{i + 1}",
                            raw={"sample": sample_name},
                        )
                    )
    return calls


def parse_longtr(vcf_path: Path) -> list[STRCall]:
    """
    Lit le VCF produit par `LongTR` (--tr-vcf, bgzippé).

    Champs utilisés (confirmés dans le README officiel gymrek-lab/LongTR) :
    INFO/MOTIF (motif du locus), INFO/END, FORMAT/GB (différence en bp de
    chaque allèle par rapport à la référence -- PAS une longueur absolue,
    contrairement à TRGT/AL). Comme pour tandem-genotypes, cette valeur
    reste cohérente en interne (comparaison max contrôle / diff patient
    au sein du même outil) même si elle n'est pas directement comparable
    aux longueurs absolues des autres outils.
    """
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        return []
    calls = []
    with pysam.VariantFile(str(vcf_path)) as vf:
        for record in vf:
            motif_field = record.info.get("MOTIF")
            if motif_field is None:
                continue
            motif = str(_first(motif_field)).split(",")[0]
            end = record.info.get("END", record.stop)

            for sample_name, sample_data in record.samples.items():
                gb = sample_data.get("GB")  # bp diff vs référence par allèle, ex: (0, 12)
                if gb is None:
                    continue
                for i, bp_diff in enumerate(gb):
                    if bp_diff is None:
                        continue
                    calls.append(
                        STRCall(
                            chrom=record.chrom,
                            start=record.pos,
                            end=int(end) if end else record.pos,
                            motif=motif,
                            size=float(bp_diff),
                            source=f"longtr_allele{i + 1}",
                            raw={"sample": sample_name},
                        )
                    )
    return calls


def _split_two_alleles(values: list[float]) -> tuple[float, float]:
    """
    Sépare une liste de longueurs par read (bp) en 2 groupes (allèle court /
    long) au niveau du plus grand écart entre valeurs triées, puis retourne
    la médiane de chaque groupe. Robuste au bruit read-level, contrairement
    à un simple min/max.
    """
    values = sorted(values)
    if len(values) == 1:
        return values[0], values[0]
    gap_idx = max(range(len(values) - 1), key=lambda i: values[i + 1] - values[i])
    low, high = values[: gap_idx + 1], values[gap_idx + 1 :]
    return median(low), median(high)


def _merge_overlapping_tg_rows(rows: list[dict]) -> list[dict]:
    """
    repeats.trf.bed liste souvent plusieurs motifs candidats (périodes
    différentes) qui se chevauchent pour le même locus TRF -- tandem-genotypes
    rapporte alors une ligne par candidat. On regroupe les lignes dont les
    intervalles se chevauchent (par chromosome) et on ne garde que celle
    couverte par le plus de reads comme représentante du locus (critère
    sans ambiguïté, contrairement à la colonne 5 -- voir parse_tandem_genotypes).
    """
    by_chrom: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_chrom[r["chrom"]].append(r)

    kept = []
    for chrom_rows in by_chrom.values():
        chrom_rows.sort(key=lambda r: r["start"])
        current_group = []
        current_end = None
        for r in chrom_rows:
            if current_group and r["start"] <= current_end:
                current_group.append(r)
                current_end = max(current_end, r["end"])
            else:
                if current_group:
                    kept.append(max(current_group, key=lambda g: len(g["values"])))
                current_group = [r]
                current_end = r["end"]
        if current_group:
            kept.append(max(current_group, key=lambda g: len(g["values"])))
    return kept


def parse_tandem_genotypes(tsv_path: Path) -> list[STRCall]:
    """
    Lit le TSV produit par `tandem-genotypes repeats.bed alignments.maf`.

    Format à 8 colonnes tab-séparées, colonne 5 (index 4, "gene name" /
    "score" selon la doc officielle -- ambigu sur nos fichiers, valeurs
    type 2.8/18.8 qui ne ressemblent pas à un nom de gène) IGNORÉE, on ne
    l'exploite pas :
      0 chrom, 1 start, 2 end, 3 motif, 4 (ignorée), 5 '.',
      6 longueurs par read (bp, séparées par virgule), 7 '.'

    La colonne 6 est une liste de longueurs de répétition mesurées par read
    individuel -- on la sépare en 2 groupes (allèle court/long, cf.
    _split_two_alleles) dont on prend la médiane comme taille d'allèle en bp
    (comparable à TRGT/AL). Les motifs candidats qui se chevauchent (issus
    de repeats.trf.bed) sont dédoublonnés en gardant celui couvert par le
    plus de reads (cf. _merge_overlapping_tg_rows).
    """
    tsv_path = Path(tsv_path)
    if not tsv_path.exists():
        return []

    rows = []
    with open(tsv_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 7:
                continue
            try:
                chrom, start, end, motif = cols[0], int(cols[1]), int(cols[2]), cols[3]
                values = [float(x) for x in cols[6].split(",") if x]
            except ValueError:
                continue
            if not values:
                continue
            rows.append({"chrom": chrom, "start": start, "end": end, "motif": motif, "values": values})

    calls = []
    for row in _merge_overlapping_tg_rows(rows):
        allele1, allele2 = _split_two_alleles(row["values"])
        for i, size in enumerate((allele1, allele2)):
            calls.append(
                STRCall(
                    chrom=row["chrom"],
                    start=row["start"],
                    end=row["end"],
                    motif=row["motif"],
                    size=size,
                    source=f"tandem_genotypes_allele{i + 1}",
                    raw={"n_reads": len(row["values"])},
                )
            )
    return calls


# ---------------------------------------------------------------------
# Matching flou par intervalle (+ marge) et motif canonique
# ---------------------------------------------------------------------

def cluster_calls(calls: list[STRCall], window: int = 25) -> list[list[STRCall]]:
    """
    Regroupe les calls en loci : deux calls sont dans le même cluster si
    leurs intervalles [start, end] se chevauchent (avec tolérance `window`
    bp de part et d'autre) ET que leurs motifs canoniques sont identiques.

    Implémenté en un seul passage trié par position (balayage), avec une
    liste de clusters "ouverts" (encore atteignables compte tenu de la
    tolérance) — pas de comparaison O(n²) sur l'ensemble du chromosome.
    """
    by_chrom: dict[str, list[STRCall]] = defaultdict(list)
    for c in calls:
        by_chrom[c.chrom].append(c)

    all_clusters: list[list[STRCall]] = []

    for chrom_calls in by_chrom.values():
        chrom_calls.sort(key=lambda c: c.start)
        open_clusters: list[dict] = []  # [{"end": int, "motifs": set, "calls": [...]}]

        for call in chrom_calls:
            cm = canonical_motif(call.motif)

            still_open = []
            for oc in open_clusters:
                if oc["end"] + window >= call.start:
                    still_open.append(oc)
                else:
                    all_clusters.append(oc["calls"])
            open_clusters = still_open

            match = next((oc for oc in open_clusters if cm in oc["motifs"]), None)
            if match is not None:
                match["calls"].append(call)
                match["end"] = max(match["end"], call.end)
                match["motifs"].add(cm)
            else:
                open_clusters.append({"end": call.end, "motifs": {cm}, "calls": [call]})

        all_clusters.extend(oc["calls"] for oc in open_clusters)

    return all_clusters


def build_locus_record(cluster: list[STRCall]) -> dict:
    chrom = cluster[0].chrom
    start = min(c.start for c in cluster)
    end = max(c.end for c in cluster)

    motif_counts = Counter(canonical_motif(c.motif) for c in cluster)
    motif = motif_counts.most_common(1)[0][0]

    sizes_by_source = {c.source: c.size for c in cluster}

    return {
        "chrom": chrom,
        "pos": start,
        "end": end,
        "motif": motif,
        "sources": sorted(sizes_by_source),
        "sizes_by_source": sizes_by_source,
    }


# ---------------------------------------------------------------------
# Écriture du VCF fusionné
# ---------------------------------------------------------------------

_VCF_HEADER = """##fileformat=VCFv4.2
##source=str-expansion-toolkit merge_to_vcf
##INFO=<ID=MOTIF,Number=1,Type=String,Description="Motif canonique du locus">
##INFO=<ID=END,Number=1,Type=Integer,Description="Fin du locus (max des sources)">
##INFO=<ID=SOURCES,Number=.,Type=String,Description="Outils/haplotypes/allèles ayant rapporté ce locus">
##INFO=<ID=SIZES,Number=.,Type=String,Description="Paires source:taille (unités différentes selon l'outil, voir README)">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def write_merged_vcf(records: list[dict], out_path: Path) -> Path:
    records = sorted(records, key=lambda r: (r["chrom"], r["pos"]))
    with open(out_path, "w") as fh:
        fh.write(_VCF_HEADER)
        for r in records:
            sizes_str = ",".join(f"{src}:{val:g}" for src, val in r["sizes_by_source"].items())
            sources_str = ",".join(r["sources"])
            info = f"MOTIF={r['motif']};END={r['end']};SOURCES={sources_str};SIZES={sizes_str}"
            fh.write(f"{r['chrom']}\t{r['pos']}\t.\tN\t.\t.\t.\t{info}\n")
    return out_path


# ---------------------------------------------------------------------
# Point d'entrée appelé par detect.merge_to_vcf
# ---------------------------------------------------------------------

def merge_tool_outputs(sample_id: str, tool_outputs: dict[str, object], outdir: Path, window: int = 25) -> Path:
    calls: list[STRCall] = []

    if "vamos" in tool_outputs:
        calls += parse_vamos(tool_outputs["vamos"])
    if "trgt" in tool_outputs:
        calls += parse_trgt(tool_outputs["trgt"])
    if "tandem-genotypes" in tool_outputs:
        calls += parse_tandem_genotypes(tool_outputs["tandem-genotypes"])
    if "longtr" in tool_outputs:
        calls += parse_longtr(tool_outputs["longtr"])

    clusters = cluster_calls(calls, window=window)
    records = [build_locus_record(c) for c in clusters]

    out_path = Path(outdir) / f"{sample_id}.merged.vcf"
    return write_merged_vcf(records, out_path)


# ---------------------------------------------------------------------
# Lecture du VCF fusionné (utilisé par build-controls / compare)
# ---------------------------------------------------------------------

# Correspondance source (ex: "vamos_hap1", "tandem_genotypes_allele2") -> outil
_TOOL_FAMILY_NAMES = {"vamos": "vamos", "trgt": "trgt", "tandem_genotypes": "tandem-genotypes", "longtr": "longtr"}


def tool_family(source: str) -> str:
    """'vamos_hap1' -> 'vamos' ; 'trgt_allele1' -> 'trgt' ; 'tandem_genotypes_allele1' -> 'tandem-genotypes'."""
    prefix = source.rsplit("_", 1)[0]
    return _TOOL_FAMILY_NAMES.get(prefix, prefix)


def parse_merged_vcf(path: Path):
    """
    Relit un VCF produit par write_merged_vcf. Génère des dicts
    {chrom, pos, end, motif, sizes_by_source} par locus.
    """
    path = Path(path)
    if not path.exists():
        return
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, _id, _ref, _alt, _qual, _filt, info = fields[:8]

            info_dict = {}
            for kv in info.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info_dict[k] = v

            sizes_by_source = {}
            for pair in info_dict.get("SIZES", "").split(","):
                if ":" not in pair:
                    continue
                src, val = pair.rsplit(":", 1)
                try:
                    sizes_by_source[src] = float(val)
                except ValueError:
                    continue

            yield {
                "chrom": chrom,
                "pos": int(pos),
                "end": int(info_dict.get("END", pos)),
                "motif": info_dict.get("MOTIF", ""),
                "sizes_by_source": sizes_by_source,
            }
