"""
`detect` subcommand.

Runs VAMOS, tandem-genotypes, and LongTR (TRGT opt-in) for one or more
samples, then merges the outputs into a single VCF per sample (see
str_toolkit/merge.py for the merge logic: fuzzy interval + canonical motif
matching, since the tools do not anchor their coordinates the same way).

Inputs: a coordinate-sorted BAM, a FASTQ, or both. A BAM alone is enough,
provided it was aligned to the reference given in config.yaml: it is passed
straight to the alignment-based tools, and the reads tandem-genotypes needs
are extracted from it. A FASTQ alone is also enough, since the alignment is
then produced once and shared. Pass realign=True to ignore an input BAM and
realign its reads instead, which a BAM from another reference build or
another aligner requires.

Prerequisite: micromamba must be installed, and the environments referenced
in config.yaml (clair3, whatshap-env, vamos, trgt, longtr, last_env,
tandem-env) must already exist on the execution machine.
"""

from __future__ import annotations

import glob
import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import NamedTuple

from str_toolkit.config import Config
from str_toolkit.utils import (
    ensure_outdir,
    read_samples_list,
    run_cmd,
    run_in_env,
    set_env_runner,
    wrap_pipeline,
)
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


def _run_maybe_env(env: str, args: list[str], *, stdout_path: Path | str | None = None,
                   shell_pipeline: str | None = None):
    """
    Runs a command inside `env` when that name is non-empty, and on the
    current PATH otherwise. Used by the steps that belong to no single tool,
    where the binary is often a system one rather than one installed in a
    tool environment.
    """
    if env:
        return run_in_env(env, args, stdout_path=stdout_path,
                          shell_pipeline=shell_pipeline)
    if shell_pipeline is not None:
        logger.info("bash -c: %s", shell_pipeline)
        return subprocess.run(["bash", "-c", wrap_pipeline(shell_pipeline)], check=True)
    if stdout_path is not None:
        logger.info("%s > %s", " ".join(args), stdout_path)
        with open(stdout_path, "w") as out_fh:
            return subprocess.run(args, stdout=out_fh, check=True)
    return run_cmd(args)


def _samtools(cfg: Config, args: list[str], *, stdout_path: Path | str | None = None,
              shell_pipeline: str | None = None):
    """
    Runs samtools either inside the environment named by `samtools_env` in
    config.yaml, or directly on the current PATH when that key is absent.
    Used for the input inspection and conversion steps, which do not belong
    to any single tool's environment.
    """
    return _run_maybe_env(cfg.samtools_env, args, stdout_path=stdout_path,
                          shell_pipeline=shell_pipeline)


def _bam_is_indexed(bam: str | Path) -> bool:
    bam = str(bam)
    candidates = [f"{bam}.bai", f"{bam}.csi", f"{bam}.crai",
                  str(Path(bam).with_suffix(".bai"))]
    return any(os.path.exists(c) for c in candidates)


def _ensure_indexed_bam(bam: str, cfg: Config, outdir: Path, threads: int, sid: str) -> Path:
    """
    Returns an indexed BAM path for an already-aligned input BAM.

    The input BAM often lives in a read-only or shared directory, so it is
    never indexed in place: when no index is found next to it, a symlink is
    created in the sample output directory and that symlink is indexed
    instead. The BAM data itself is never copied.
    """
    if _bam_is_indexed(bam):
        return Path(bam)

    link = outdir / f"{sid}.input.bam"
    if not link.exists():
        logger.info("[%s] input BAM has no index, symlinking into %s", sid, outdir)
        link.symlink_to(os.path.abspath(bam))
    if not _bam_is_indexed(link):
        logger.info("[%s] indexing %s", sid, link.name)
        _samtools(cfg, ["samtools", "index", "-@", str(threads), str(link)])
    return link


# ---------------------------------------------------------------------
# VAMOS: clair3 (phasing) -> whatshap haplotag/split -> vamos --contig x2
# ---------------------------------------------------------------------

def run_vamos(sample: Sample, cfg: Config, outdir: Path, threads: int,
              realign: bool = False) -> dict[str, Path]:
    sid = sample.sample_id
    cvamos = cfg.vamos

    # Both final VCFs already present: nothing to do, and in particular no
    # need for the intermediates. The per-step skips below are sequential, so
    # without this short circuit a deleted haplotagged BAM would be rebuilt
    # even though the VCFs it leads to are already there. Those intermediates
    # are large (a haplotagged BAM is about the size of the input), and
    # deleting them after a finished run has to stay safe.
    hap_vcfs_done = {
        hap: outdir / f"{sid}_assembly.{hap}.vcf" for hap in ("hap1", "hap2")
    }
    if all(v.exists() for v in hap_vcfs_done.values()):
        logger.info("[%s] VAMOS: both haplotype VCFs already present, skipping", sid)
        return hap_vcfs_done

    # clair3 needs an indexed, aligned BAM. When only a FASTQ is available,
    # or when the caller asked for a realignment, fall back to the shared
    # minimap2 step rather than failing.
    if sample.bam_path and not realign:
        bam = str(_ensure_indexed_bam(sample.bam_path, cfg, outdir, threads, sid))
    else:
        bam = str(_ensure_ont_sorted_bam(
            sample, cfg.longtr.mmi or cfg.trgt.mmi, cfg.longtr.env or cfg.trgt.env,
            outdir, threads, "VAMOS", cfg, realign,
        ))
        if not _bam_is_indexed(bam):
            _samtools(cfg, ["samtools", "index", "-@", str(threads), bam])

    phased_vcf = outdir / "phased_merge_output.vcf.gz"
    if not phased_vcf.exists():
        logger.info("[%s] VAMOS: clair3 (phasing)", sid)
        _run_maybe_env(
            cvamos.env_clair3,
            [
                cvamos.bin_clair3,
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
        _run_maybe_env(
            cvamos.env_whatshap,
            [
                cvamos.bin_whatshap, "haplotag",
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
        _run_maybe_env(
            cvamos.env_whatshap,
            [
                cvamos.bin_whatshap, "split",
                "--output-h1", str(h1),
                "--output-h2", str(h2),
                bam,
                str(haplotype_tsv),
            ],
        )
        _samtools(cfg, ["samtools", "index", str(h1)])
        _samtools(cfg, ["samtools", "index", str(h2)])
    else:
        logger.info("[%s] VAMOS: h1/h2 bam already present, skipping", sid)

    hap_vcfs = {}
    for hap, hap_bam in (("hap1", h1), ("hap2", h2)):
        hap_vcf = outdir / f"{sid}_assembly.{hap}.vcf"
        if not hap_vcf.exists():
            logger.info("[%s] VAMOS: vamos --contig (%s)", sid, hap)
            _run_maybe_env(
                cvamos.env_vamos,
                [
                    cvamos.bin_vamos, "--contig",
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
# Shared inputs (an aligned BAM, and a FASTQ), reused across tools.
#
# TRGT and LongTR both read an already-aligned BAM, so the alignment is
# done at most once per job. An input BAM supplied by the caller is reused
# as is rather than realigned, provided it was aligned to the same
# reference as `reference` in config.yaml. Use realign=True to force a
# fresh minimap2 alignment from the FASTQ instead, which is what an input
# BAM aligned to another build, or with an unsuitable preset, requires.
# ---------------------------------------------------------------------

def _ensure_ont_sorted_bam(
    sample: Sample,
    mmi: str,
    align_env: str,
    outdir: Path,
    threads: int,
    tool_label: str,
    cfg: Config,
    realign: bool = False,
) -> Path:
    sid = sample.sample_id
    existing = sorted(glob.glob(str(outdir / f"{sid}*.sorted.bam")))
    if existing:
        return Path(existing[0])

    if sample.bam_path and not realign:
        logger.info("[%s] %s: reusing the input BAM (no realignment)", sid, tool_label)
        return _ensure_indexed_bam(sample.bam_path, cfg, outdir, threads, sid)

    fastq = _require(sample.fastq_path, sid, "fastq_path", tool_label)
    sorted_bam = outdir / f"{sid}.sorted.bam"
    # align_env in config.yaml wins, so minimap2 does not have to be installed
    # inside the genotyping tool's own environment. Empty means the PATH.
    env = cfg.align_env if cfg.align_env else align_env
    logger.info("[%s] %s: minimap2 (map-ont) align + sort, env=%s",
                sid, tool_label, env or "PATH")
    tmp_bam = outdir / f"{sid}.sorted.bam.tmp"
    _run_maybe_env(
        env,  # must provide minimap2 + samtools
        [],
        shell_pipeline=(
            f"minimap2 -t {threads} -ax map-ont -Y {mmi} {fastq} "
            f"| samtools sort -@ {threads} -o {tmp_bam}"
        ),
    )
    os.replace(tmp_bam, sorted_bam)
    return sorted_bam


def _ensure_fastq(sample: Sample, cfg: Config, outdir: Path, threads: int, tool_label: str) -> str:
    """
    Returns a FASTQ for the sample, extracting one from the input BAM when
    none was supplied. tandem-genotypes needs reads rather than an
    alignment, since it realigns them with LAST.

    Extraction keeps the MM/ML methylation tags and, for a BAM that is
    already aligned, drops secondary and supplementary records so that each
    read appears exactly once.
    """
    if sample.fastq_path:
        return sample.fastq_path

    sid = sample.sample_id
    bam = _require(sample.bam_path, sid, "fastq_path or bam_path", tool_label)

    fastq_out = outdir / f"{sid}.reads.fastq.gz"
    if fastq_out.exists():
        logger.info("[%s] %s: %s already present, skipping extraction", sid, tool_label, fastq_out.name)
        return str(fastq_out)

    logger.info("[%s] %s: no FASTQ supplied, extracting reads from the input BAM", sid, tool_label)
    tmp_out = outdir / f"{sid}.reads.fastq.gz.tmp"
    # Whole-genome extraction is compression-bound, and single-threaded gzip
    # turns a step of minutes into one of hours, so pigz is used when present.
    compress = f"pigz -p {threads} -c" if shutil.which("pigz") else "gzip -c"
    _samtools(
        cfg,
        [],
        shell_pipeline=(
            f"samtools fastq -@ {threads} -T MM,ML -F 0x900 {bam} "
            f"| {compress} > {tmp_out}"
        ),
    )
    os.replace(tmp_out, fastq_out)
    return str(fastq_out)


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

def run_trgt(sample: Sample, cfg: Config, outdir: Path, threads: int,
             realign: bool = False) -> Path:
    sid = sample.sample_id
    ctrgt = cfg.trgt

    sorted_bam = _ensure_ont_sorted_bam(
        sample, ctrgt.mmi, ctrgt.env, outdir, threads, "TRGT", cfg, realign
    )

    if not _bam_is_indexed(sorted_bam):
        logger.info("[%s] TRGT: indexing bam", sid)
        _samtools(cfg, ["samtools", "index", "-@", str(threads), str(sorted_bam)])

    out_prefix = outdir / f"{sid}.trgt"
    out_vcf = Path(f"{out_prefix}.vcf.gz")
    if not out_vcf.exists():
        logger.info("[%s] TRGT: trgt genotype", sid)
        _run_maybe_env(
            ctrgt.env,
            [
                ctrgt.bin, "genotype",
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

def run_longtr(sample: Sample, cfg: Config, outdir: Path, threads: int,
               realign: bool = False) -> Path:
    sid = sample.sample_id
    clongtr = cfg.longtr

    sorted_bam = _ensure_ont_sorted_bam(
        sample, clongtr.mmi, clongtr.env, outdir, threads, "LongTR", cfg, realign
    )

    if not _bam_is_indexed(sorted_bam):
        logger.info("[%s] LongTR: indexing bam", sid)
        _samtools(cfg, ["samtools", "index", "-@", str(threads), str(sorted_bam)])

    out_vcf = outdir / f"{sid}.longtr.vcf.gz"
    if not out_vcf.exists():
        logger.info("[%s] LongTR: genotyping", sid)
        # LongTR has no native multi-threading; --bam-samps/--bam-libs avoids
        # relying on correct @RG tags in the BAM produced by minimap2.
        _run_maybe_env(
            clongtr.env,
            [
                clongtr.bin,
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

def run_tandem_genotypes(sample: Sample, cfg: Config, outdir: Path, threads: int,
                         realign: bool = False) -> Path:
    sid = sample.sample_id
    ctg = cfg.tandem_genotypes

    tsv_out = outdir / f"{sid}.tandem_genotypes.tsv"
    if tsv_out.exists():
        logger.info("[%s] tandem-genotypes: already done, skipping", sid)
        return tsv_out

    fastq = _ensure_fastq(sample, cfg, outdir, threads, "tandem-genotypes")

    # LAST memory grows with the thread count and with the query batch size,
    # and lastal on a whole ONT genome has exhausted RAM and swap at the
    # defaults. `last_threads` and `lastal_batch_size` bound both; see
    # config.example.yaml.
    last_threads = ctg.last_threads if ctg.last_threads > 0 else threads

    par_file = outdir / f"{sid}_reads.par"
    if par_file.exists():
        logger.info("[%s] tandem-genotypes: last-train parameters already present", sid)
    else:
        logger.info("[%s] tandem-genotypes: last-train (%d threads)", sid, last_threads)
        tmp_par = outdir / f"{sid}_reads.par.tmp"
        _run_maybe_env(
            ctg.env_last,
            [ctg.bin_last_train, "-P", str(last_threads), "-Q0", ctg.last_ref_db, fastq],
            stdout_path=tmp_par,
        )
        os.replace(tmp_par, par_file)

    maf_file = outdir / f"{sid}_alignments.maf"
    if maf_file.exists():
        logger.info("[%s] tandem-genotypes: alignments already present", sid)
    else:
        logger.info("[%s] tandem-genotypes: lastal | last-split (%d threads, batch %s)",
                    sid, last_threads, ctg.lastal_batch_size or "default")
        # Written under a temporary name and renamed only on a clean exit, so a
        # killed lastal can never leave a truncated MAF that the check above
        # would then accept as finished.
        tmp_maf = outdir / f"{sid}_alignments.maf.tmp"
        batch = f"-i {ctg.lastal_batch_size} " if ctg.lastal_batch_size else ""
        extra = f"{ctg.lastal_extra_args} " if ctg.lastal_extra_args else ""
        _run_maybe_env(
            ctg.env_last,
            [],
            shell_pipeline=(
                f"{ctg.bin_lastal} -P{last_threads} {batch}{extra}--split "
                f"-p {par_file} {ctg.last_ref_db} {fastq} "
                f"| {ctg.bin_last_split} -m1 > {tmp_maf}"
            ),
        )
        os.replace(tmp_maf, maf_file)

    logger.info("[%s] tandem-genotypes: genotyping", sid)
    tmp_tsv = outdir / f"{sid}.tandem_genotypes.tsv.tmp"
    _run_maybe_env(
        ctg.env_tandem,
        [ctg.bin_tandem_genotypes, ctg.repeats_bed, str(maf_file)],
        stdout_path=tmp_tsv,
    )
    os.replace(tmp_tsv, tsv_out)

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
    sample: Sample, cfg: Config, outdir: Path, tools: list[str], threads: int,
    realign: bool = False,
) -> Path:
    sample_outdir = outdir / sample.sample_id
    ensure_outdir(sample_outdir)

    tool_outputs: dict[str, object] = {}
    for tool_name in tools:
        logger.info("Sample %s: running %s", sample.sample_id, tool_name)
        runner = TOOL_RUNNERS[tool_name]
        tool_outputs[tool_name] = runner(sample, cfg, sample_outdir, threads, realign)

    final_vcf = merge_to_vcf(sample, tool_outputs, sample_outdir)
    logger.info("Sample %s: final VCF -> %s", sample.sample_id, final_vcf)
    return final_vcf


def run(args) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    cfg = Config.from_yaml(args.config)
    set_env_runner(cfg.env_runner or None)
    outdir = Path(args.outdir)
    ensure_outdir(outdir)

    if args.sample:
        if not args.bam and not args.fastq:
            raise SystemExit("--bam and/or --fastq is required when --sample is used")
        samples = [Sample(sample_id=args.sample, bam_path=args.bam, fastq_path=args.fastq)]
    else:
        samples = [Sample(**row) for row in read_samples_list(args.samples_list)]

    produced_vcfs = []
    realign = getattr(args, "realign", False)
    for sample in samples:
        vcf_path = detect_one_sample(sample, cfg, outdir, args.tools, args.threads, realign)
        produced_vcfs.append(vcf_path)

    logger.info("Done: %d VCFs produced in %s", len(produced_vcfs), outdir)
    return 0
