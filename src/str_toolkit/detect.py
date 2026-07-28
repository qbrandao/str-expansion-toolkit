"""
Sous-commande `detect`.

Lance VAMOS, TRGT et tandem-genotypes pour un ou plusieurs patients,
puis fusionne les 3 sorties en un unique VCF par patient (voir
str_toolkit/merge.py pour la logique de fusion : matching flou par
intervalle + motif canonique, car les 3 outils n'ancrent pas leurs
coordonnées de la même façon).

Prérequis : micromamba doit être installé et les environnements référencés
dans le config.yaml (clair3, whatshap-env, vamos, trgt, last_env, tandem-env)
doivent déjà exister sur la machine d'exécution.
"""

from __future__ import annotations

import glob
import logging
from pathlib import Path
from typing import NamedTuple

from str_toolkit.config import Config
from str_toolkit.utils import ensure_outdir, read_samples_list, run_in_env
from str_toolkit import merge

logger = logging.getLogger(__name__)


class Sample(NamedTuple):
    sample_id: str
    bam_path: str | None = None    # BAM déjà aligné (utilisé par VAMOS/clair3)
    fastq_path: str | None = None  # fastq(.gz) brut fusionné (utilisé par TRGT et tandem-genotypes)


def _require(value, sample_id: str, field_name: str, tool: str) -> str:
    if not value:
        raise SystemExit(
            f"Sample {sample_id}: '{field_name}' est requis pour lancer {tool} "
            f"(à fournir via --bam/--fastq ou la colonne correspondante du TSV)."
        )
    return value


# ---------------------------------------------------------------------
# VAMOS : clair3 (phasing) -> whatshap haplotag/split -> vamos --contig x2
# ---------------------------------------------------------------------

def run_vamos(sample: Sample, cfg: Config, outdir: Path, threads: int) -> dict[str, Path]:
    bam = _require(sample.bam_path, sample.sample_id, "bam_path", "VAMOS")
    sid = sample.sample_id
    cvamos = cfg.vamos

    phased_vcf = outdir / "phased_merge_output.vcf.gz"
    if not phased_vcf.exists():
        logger.info("[%s] VAMOS: clair3 (phasing)", sid)
        run_in_env(
            cvamos.env_clair3,
            [
                "run_clair3.sh",
                f"--bam_fn={bam}",
                f"--ref_fn={cfg.reference}",
                f"--threads={threads}",
                "--platform=ont",
                f"--model_path={cvamos.model_prefix}",
                f"--output={outdir}",
                "--use_whatshap_for_final_output_phasing",
            ],
        )
    else:
        logger.info("[%s] VAMOS: clair3 déjà fait, skip", sid)

    haplotagged_bam = outdir / f"{sid}_haplotagged.bam"
    haplotype_tsv = outdir / f"{sid}_haplotype.tsv"
    if not haplotagged_bam.exists():
        logger.info("[%s] VAMOS: whatshap haplotag", sid)
        run_in_env(
            cvamos.env_whatshap,
            [
                "whatshap", "haplotag",
                "-o", str(haplotagged_bam),
                "--reference", cfg.reference,
                str(phased_vcf),
                bam,
                "--output-haplotag-list", str(haplotype_tsv),
                "--ignore-read-groups",
                f"--output-threads={threads}",
            ],
        )
    else:
        logger.info("[%s] VAMOS: haplotagged bam déjà présent, skip", sid)

    h1 = outdir / f"{sid}_h1.bam"
    h2 = outdir / f"{sid}_h2.bam"
    if not (h1.exists() and h2.exists()):
        logger.info("[%s] VAMOS: whatshap split", sid)
        run_in_env(
            cvamos.env_whatshap,
            [
                "whatshap", "split",
                "--output-h1", str(h1),
                "--output-h2", str(h2),
                bam,
                str(haplotype_tsv),
            ],
        )
        run_in_env(cvamos.env_whatshap, ["samtools", "index", str(h1)])
        run_in_env(cvamos.env_whatshap, ["samtools", "index", str(h2)])
    else:
        logger.info("[%s] VAMOS: h1/h2 bam déjà présents, skip", sid)

    hap_vcfs = {}
    for hap, hap_bam in (("hap1", h1), ("hap2", h2)):
        hap_vcf = outdir / f"{sid}_assembly.{hap}.vcf"
        if not hap_vcf.exists():
            logger.info("[%s] VAMOS: vamos --contig (%s)", sid, hap)
            run_in_env(
                cvamos.env_vamos,
                [
                    "vamos", "--contig",
                    "-b", str(hap_bam),
                    "-r", cvamos.catalog,
                    "-s", sid,
                    "-o", str(hap_vcf),
                    "-t", str(threads),
                ],
            )
        else:
            logger.info("[%s] VAMOS: %s déjà présent, skip", sid, hap_vcf.name)
        hap_vcfs[hap] = hap_vcf

    return hap_vcfs


# ---------------------------------------------------------------------
# Alignement partagé (minimap2 ONT + tri + indexation), réutilisé par TRGT
# et LongTR : les deux outils lisent un BAM/CRAM déjà aligné, inutile de
# ré-aligner deux fois si les deux tournent dans le même run.
# ---------------------------------------------------------------------

def _ensure_ont_sorted_bam(sample: Sample, mmi: str, align_env: str, outdir: Path, threads: int, tool_label: str) -> Path:
    sid = sample.sample_id
    existing = sorted(glob.glob(str(outdir / f"{sid}*.sorted.bam")))
    if existing:
        return Path(existing[0])

    fastq = _require(sample.fastq_path, sid, "fastq_path", tool_label)
    sorted_bam = outdir / f"{sid}.sorted.bam"
    logger.info("[%s] %s: minimap2 (map-ont) align + sort", sid, tool_label)
    run_in_env(
        align_env,  # doit fournir minimap2 + samtools
        [],
        shell_pipeline=(
            f"minimap2 -t {threads} -ax map-ont -Y {mmi} {fastq} "
            f"| samtools sort -@ {threads} -o {sorted_bam}"
        ),
    )
    return sorted_bam


# ---------------------------------------------------------------------
# TRGT : minimap2 (align + sort) -> trgt genotype
#
# ATTENTION : TRGT est conçu pour des reads PacBio HiFi et n'a pas de
# support officiel pour les données ONT (cf. Aliyev et al. 2026, bioRxiv,
# qui exclut explicitement TRGT des benchmarks ONT pour cette raison). Ce
# n'est donc PAS un outil par défaut de `detect` -- il ne tourne que si
# explicitement demandé via `--tools ... trgt ...`. À documenter comme
# usage hors cadre officiel dans toute publication utilisant ces résultats.
# ---------------------------------------------------------------------

def run_trgt(sample: Sample, cfg: Config, outdir: Path, threads: int) -> Path:
    sid = sample.sample_id
    ctrgt = cfg.trgt

    sorted_bam = _ensure_ont_sorted_bam(sample, ctrgt.mmi, ctrgt.env, outdir, threads, "TRGT")

    bai = Path(f"{sorted_bam}.bai")
    if not bai.exists():
        logger.info("[%s] TRGT: indexation bam", sid)
        run_in_env(ctrgt.env, ["samtools", "index", "-@", str(threads), str(sorted_bam)])

    out_prefix = outdir / f"{sid}.trgt"
    out_vcf = Path(f"{out_prefix}.vcf.gz")
    if not out_vcf.exists():
        logger.info("[%s] TRGT: trgt genotype", sid)
        run_in_env(
            ctrgt.env,
            [
                "trgt", "genotype",
                "--threads", str(threads),
                "--reads", str(sorted_bam),
                "--genome", cfg.reference,
                "--repeats", ctrgt.repeats_bed,
                "--output-prefix", str(out_prefix),
            ],
        )
    else:
        logger.info("[%s] TRGT: %s déjà présent, skip", sid, out_vcf.name)

    return out_vcf


# ---------------------------------------------------------------------
# LongTR : minimap2 (align + sort, partagé avec TRGT) -> LongTR
#
# Outil ONT-natif (HipSTR adapté aux reads longs, PacBio HiFi ET ONT).
# Choisi comme 3e outil par défaut à la place de TRGT pour les données ONT
# -- meilleure concordance avec les assemblages, mais nécessite des reads
# de bonne qualité/profondeur suffisante (--min-reads=10 par défaut).
# ---------------------------------------------------------------------

def run_longtr(sample: Sample, cfg: Config, outdir: Path, threads: int) -> Path:
    sid = sample.sample_id
    clongtr = cfg.longtr

    sorted_bam = _ensure_ont_sorted_bam(sample, clongtr.mmi, clongtr.env, outdir, threads, "LongTR")

    bai = Path(f"{sorted_bam}.bai")
    if not bai.exists():
        logger.info("[%s] LongTR: indexation bam", sid)
        run_in_env(clongtr.env, ["samtools", "index", "-@", str(threads), str(sorted_bam)])

    out_vcf = outdir / f"{sid}.longtr.vcf.gz"
    if not out_vcf.exists():
        logger.info("[%s] LongTR: genotyping", sid)
        # LongTR n'a pas de multi-threading natif ; --bam-samps/--bam-libs
        # évite de dépendre de tags @RG corrects dans le BAM produit par minimap2.
        run_in_env(
            clongtr.env,
            [
                "LongTR",
                "--bams", str(sorted_bam),
                "--fasta", cfg.reference,
                "--regions", clongtr.regions_bed,
                "--tr-vcf", str(out_vcf),
                "--bam-samps", sid,
                "--bam-libs", f"{sid}_lib",
            ],
        )
    else:
        logger.info("[%s] LongTR: %s déjà présent, skip", sid, out_vcf.name)

    return out_vcf


# ---------------------------------------------------------------------
# tandem-genotypes : last-train -> lastal | last-split -> tandem-genotypes
# ---------------------------------------------------------------------

def run_tandem_genotypes(sample: Sample, cfg: Config, outdir: Path, threads: int) -> Path:
    sid = sample.sample_id
    ctg = cfg.tandem_genotypes

    tsv_out = outdir / f"{sid}.tandem_genotypes.tsv"
    if tsv_out.exists():
        logger.info("[%s] tandem-genotypes: déjà fait, skip", sid)
        return tsv_out

    fastq = _require(sample.fastq_path, sid, "fastq_path", "tandem-genotypes")

    par_file = outdir / f"{sid}_reads.par"
    logger.info("[%s] tandem-genotypes: last-train", sid)
    run_in_env(
        ctg.env_last,
        ["last-train", "-P", str(threads), "-Q0", ctg.last_ref_db, fastq],
        stdout_path=par_file,
    )

    maf_file = outdir / f"{sid}_alignments.maf"
    logger.info("[%s] tandem-genotypes: lastal | last-split", sid)
    run_in_env(
        ctg.env_last,
        [],
        shell_pipeline=(
            f"lastal -P{threads} --split -p {par_file} {ctg.last_ref_db} {fastq} "
            f"| last-split -m1 > {maf_file}"
        ),
    )

    logger.info("[%s] tandem-genotypes: genotyping", sid)
    run_in_env(
        ctg.env_tandem,
        ["tandem-genotypes", ctg.repeats_bed, str(maf_file)],
        stdout_path=tsv_out,
    )

    return tsv_out


TOOL_RUNNERS = {
    "vamos": run_vamos,
    "trgt": run_trgt,
    "tandem-genotypes": run_tandem_genotypes,
    "longtr": run_longtr,
}


# ---------------------------------------------------------------------
# Fusion des 3 sorties en un VCF unique (str_toolkit.merge)
# ---------------------------------------------------------------------

def merge_to_vcf(sample: Sample, tool_outputs: dict[str, object], outdir: Path) -> Path:
    """
    Fusionne les sorties des outils lancés en un unique VCF par patient, via
    un matching flou par intervalle + motif canonique (voir str_toolkit/merge.py
    pour le détail de l'algorithme et les limites, notamment sur les unités
    de taille qui diffèrent selon l'outil).
    """
    return merge.merge_tool_outputs(sample.sample_id, tool_outputs, outdir)


# ---------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------

def detect_one_sample(
    sample: Sample, cfg: Config, outdir: Path, tools: list[str], threads: int
) -> Path:
    sample_outdir = outdir / sample.sample_id
    ensure_outdir(sample_outdir)

    tool_outputs: dict[str, object] = {}
    for tool_name in tools:
        logger.info("Sample %s: running %s", sample.sample_id, tool_name)
        runner = TOOL_RUNNERS[tool_name]
        tool_outputs[tool_name] = runner(sample, cfg, sample_outdir, threads)

    final_vcf = merge_to_vcf(sample, tool_outputs, sample_outdir)
    logger.info("Sample %s: VCF final -> %s", sample.sample_id, final_vcf)
    return final_vcf


def run(args) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    cfg = Config.from_yaml(args.config)
    outdir = Path(args.outdir)
    ensure_outdir(outdir)

    if args.sample:
        if not args.bam and not args.fastq:
            raise SystemExit("--bam et/ou --fastq requis quand --sample est utilisé")
        samples = [Sample(sample_id=args.sample, bam_path=args.bam, fastq_path=args.fastq)]
    else:
        samples = [Sample(**row) for row in read_samples_list(args.samples_list)]

    produced_vcfs = []
    for sample in samples:
        vcf_path = detect_one_sample(sample, cfg, outdir, args.tools, args.threads)
        produced_vcfs.append(vcf_path)

    logger.info("Terminé : %d VCF produits dans %s", len(produced_vcfs), outdir)
    return 0
