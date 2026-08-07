"""
`detect` subcommand.

Runs VAMOS, tandem-genotypes, and LongTR (TRGT opt-in) for one or more
samples, then merges the outputs into a single VCF per sample (see
str_toolkit/merge.py for the merge logic: fuzzy interval + canonical motif
matching, since the tools do not anchor their coordinates the same way).

Prerequisite: micromamba must be installed, and the environments referenced
in config.yaml (clair3, whatshap-env, vamos, trgt, longtr, last_env,
tandem-env) must already exist on the execution machine.
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
    bam_path: str | None = None    # already-aligned BAM (used by VAMOS/clair3)
    fastq_path: str | None = None  # raw merged fastq(.gz) (used by TRGT/LongTR and tandem-genotypes)


def _require(value, sample_id: str, field_name: str, tool: str) -> str:
    if not value:
        raise SystemExit(
            f"Sample {sample_id}: '{field_name}' is required to run {tool} "
            f"(provide it via --bam/--fastq or the corresponding TSV column)."
        )
    return value


# ---------------------------------------------------------------------
# VAMOS: clair3 (phasing) -> whatshap haplotag/split -> vamos --contig x2
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
        logger.info("[%s] VAMOS: clair3 already done, skipping", sid)

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
        logger.info("[%s] VAMOS: haplotagged bam already present, skipping", sid)

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
        logger.info("[%s] VAMOS: h1/h2 bam already present, skipping", sid)

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
            logger.info("[%s] VAMOS: %s already present, skipping", sid, hap_vcf.name)
        hap_vcfs[hap] = hap_vcf

    return hap_vcfs


# ---------------------------------------------------------------------
# Shared alignment (ONT minimap2 + sort + index), reused by TRGT and
# LongTR: both tools read an already-aligned BAM/CRAM, so there is no need
# to align twice if both run in the same job.
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
        align_env,  # must provide minimap2 + samtools
        [],
        shell_pipeline=(
            f"minimap2 -t {threads} -ax map-ont -Y {mmi} {fastq} "
            f"| samtools sort -@ {threads} -o {sorted_bam}"
        ),
    )
    return sorted_bam


# ---------------------------------------------------------------------
# TRGT: minimap2 (align + sort) -> trgt genotype
#
# NOTE: TRGT is designed for PacBio HiFi reads and has no official support
# for ONT data (see Aliyev et al. 2026, bioRxiv, which explicitly excludes
# TRGT from ONT benchmarks for this reason). It is therefore NOT a default
# tool in `detect` -- it only runs if explicitly requested via
# `--tools ... trgt ...`. Any publication using these results should
# document this as an off-label use.
# ---------------------------------------------------------------------

def run_trgt(sample: Sample, cfg: Config, outdir: Path, threads: int) -> Path:
    sid = sample.sample_id
    ctrgt = cfg.trgt

    sorted_bam = _ensure_ont_sorted_bam(sample, ctrgt.mmi, ctrgt.env, outdir, threads, "TRGT")

    bai = Path(f"{sorted_bam}.bai")
    if not bai.exists():
        logger.info("[%s] TRGT: indexing bam", sid)
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
        logger.info("[%s] TRGT: %s already present, skipping", sid, out_vcf.name)

    return out_vcf


# ---------------------------------------------------------------------
# LongTR: minimap2 (align + sort, shared with TRGT) -> LongTR
#
# ONT-native tool (a long-read adaptation of HipSTR, supporting both
# PacBio HiFi and ONT). Chosen as the third default tool in place of TRGT
# for ONT data -- better concordance with assemblies, but requires
# sufficient read quality/depth (--min-reads=10 by default).
# ---------------------------------------------------------------------

def run_longtr(sample: Sample, cfg: Config, outdir: Path, threads: int) -> Path:
    sid = sample.sample_id
    clongtr = cfg.longtr

    sorted_bam = _ensure_ont_sorted_bam(sample, clongtr.mmi, clongtr.env, outdir, threads, "LongTR")

    bai = Path(f"{sorted_bam}.bai")
    if not bai.exists():
        logger.info("[%s] LongTR: indexing bam", sid)
        run_in_env(clongtr.env, ["samtools", "index", "-@", str(threads), str(sorted_bam)])

    out_vcf = outdir / f"{sid}.longtr.vcf.gz"
    if not out_vcf.exists():
        logger.info("[%s] LongTR: genotyping", sid)
        # LongTR has no native multi-threading; --bam-samps/--bam-libs avoids
        # relying on correct @RG tags in the BAM produced by minimap2.
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
        logger.info("[%s] LongTR: %s already present, skipping", sid, out_vcf.name)

    return out_vcf


# ---------------------------------------------------------------------
# tandem-genotypes: last-train -> lastal | last-split -> tandem-genotypes
# ---------------------------------------------------------------------

def run_tandem_genotypes(sample: Sample, cfg: Config, outdir: Path, threads: int) -> Path:
    sid = sample.sample_id
    ctg = cfg.tandem_genotypes

    tsv_out = outdir / f"{sid}.tandem_genotypes.tsv"
    if tsv_out.exists():
        logger.info("[%s] tandem-genotypes: already done, skipping", sid)
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
# Merge the tool outputs into a single VCF (str_toolkit.merge)
# ---------------------------------------------------------------------

def merge_to_vcf(sample: Sample, tool_outputs: dict[str, object], outdir: Path) -> Path:
    """
    Merges the outputs of the tools that were run into a single VCF per
    sample, via fuzzy interval + canonical motif matching (see
    str_toolkit/merge.py for the algorithm and its limitations, in
    particular the size units, which differ across tools).
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
    logger.info("Sample %s: final VCF -> %s", sample.sample_id, final_vcf)
    return final_vcf


def run(args) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    cfg = Config.from_yaml(args.config)
    outdir = Path(args.outdir)
    ensure_outdir(outdir)

    if args.sample:
        if not args.bam and not args.fastq:
            raise SystemExit("--bam and/or --fastq is required when --sample is used")
        samples = [Sample(sample_id=args.sample, bam_path=args.bam, fastq_path=args.fastq)]
    else:
        samples = [Sample(**row) for row in read_samples_list(args.samples_list)]

    produced_vcfs = []
    for sample in samples:
        vcf_path = detect_one_sample(sample, cfg, outdir, args.tools, args.threads)
        produced_vcfs.append(vcf_path)

    logger.info("Done: %d VCFs produced in %s", len(produced_vcfs), outdir)
    return 0
