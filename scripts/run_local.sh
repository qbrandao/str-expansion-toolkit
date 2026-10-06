#!/bin/bash
# Run `str-toolkit detect` on a machine without SLURM.
#
# The .sbatch scripts here rely on SLURM_CPUS_PER_TASK and SLURM_SUBMIT_DIR
# and cannot be used outside a cluster. This script does the same work with
# plain shell: it sizes the thread count from the machine, checks the
# environments and inputs before starting anything slow, and runs in the
# background under nohup so the run survives the terminal closing.
#
# Usage:
#   ./scripts/run_local.sh --sample ID --bam FILE [options]
#
# Options:
#   --sample ID        sample identifier, used for output names (required)
#   --bam FILE         aligned, coordinate-sorted BAM
#   --fastq FILE       reads; preferred over extracting them from the BAM
#   --config FILE      config.yaml (default: config.yaml)
#   --tools "A B"      tools to run (default: the package default set)
#   --regions FILE     override longtr.regions_bed, for a quick reduced run
#   --threads N        default: all cores minus one, capped at 16
#   --outdir DIR       default: results/local
#   --realign          realign from the FASTQ instead of reusing the BAM
#   --background       detach under nohup and return immediately
#   --status           report what has already been produced, then exit
#   --preflight        run the checks only, then exit
#   --dry-run          print the command that would run, then exit
#
# Examples:
#   # fast first pass: LongTR on chr9 only, in the background
#   awk '$1=="chr9"' /path/longtr.regions.bed > quick_chr9.bed
#   ./scripts/run_local.sh --sample 1312_DFT --bam /path/sample.bam \
#       --tools longtr --regions quick_chr9.bed --background
#
#   # full run, all default tools
#   ./scripts/run_local.sh --sample 1312_DFT --bam /path/sample.bam --background
#
#   # follow progress
#   ./scripts/run_local.sh --sample 1312_DFT --status
#   tail -f logs/local_1312_DFT.log

set -uo pipefail

SAMPLE=""; BAM=""; FASTQ=""; CONFIG="config.yaml"; TOOLS=""; REGIONS=""
RUNNER=""
THREADS=""; OUTDIR="results/local"; REALIGN=0
BACKGROUND=0; STATUS=0; PREFLIGHT_ONLY=0; DRYRUN=0

die() { echo "ERROR: $*" >&2; exit 1; }
note() { echo "[$(date +%H:%M:%S)] $*"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample)     SAMPLE="${2:?}"; shift 2 ;;
    --bam)        BAM="${2:?}"; shift 2 ;;
    --fastq)      FASTQ="${2:?}"; shift 2 ;;
    --config)     CONFIG="${2:?}"; shift 2 ;;
    --tools)      TOOLS="${2:?}"; shift 2 ;;
    --regions)    REGIONS="${2:?}"; shift 2 ;;
    --threads)    THREADS="${2:?}"; shift 2 ;;
    --outdir)     OUTDIR="${2:?}"; shift 2 ;;
    --realign)    REALIGN=1; shift ;;
    --background) BACKGROUND=1; shift ;;
    --status)     STATUS=1; shift ;;
    --preflight)  PREFLIGHT_ONLY=1; shift ;;
    --dry-run)    DRYRUN=1; shift ;;
    -h|--help)    sed -n '2,40p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *)            die "unknown option: $1" ;;
  esac
done

[[ -n "${SAMPLE}" ]] || die "--sample is required"

SAMPLE_OUT="${OUTDIR}/${SAMPLE}"
LOG="logs/local_${SAMPLE}.log"
PIDFILE="logs/local_${SAMPLE}.pid"

# --- status -----------------------------------------------------------
if [[ "${STATUS}" -eq 1 ]]; then
  echo "sample    : ${SAMPLE}"
  echo "output dir: ${SAMPLE_OUT}"
  if [[ -f "${PIDFILE}" ]] && kill -0 "$(cat "${PIDFILE}")" 2>/dev/null; then
    echo "state     : RUNNING (pid $(cat "${PIDFILE}"))"
  elif [[ -f "${PIDFILE}" ]]; then
    echo "state     : not running, and the pid file was left behind, so the run"
    echo "            was killed rather than finishing. Relaunch to resume."
  elif [[ -f "${LOG}" ]]; then
    echo "state     : not running. A previous run finished, see the exit code in"
    echo "            ${LOG}"
  else
    echo "state     : no run recorded"
  fi
  echo
  echo "expected outputs:"
  for f in "phased_merge_output.vcf.gz:clair3 phasing (VAMOS)" \
           "${SAMPLE}_haplotagged.bam:whatshap haplotag (VAMOS)" \
           "${SAMPLE}_h1.bam:haplotype 1 (VAMOS)" \
           "${SAMPLE}_h2.bam:haplotype 2 (VAMOS)" \
           "${SAMPLE}_assembly.hap1.vcf:vamos hap1" \
           "${SAMPLE}_assembly.hap2.vcf:vamos hap2" \
           "${SAMPLE}.reads.fastq.gz:reads extracted from the BAM" \
           "${SAMPLE}_reads.par:last-train parameters" \
           "${SAMPLE}_alignments.maf:LAST alignments" \
           "${SAMPLE}.tandem_genotypes.tsv:tandem-genotypes" \
           "${SAMPLE}.longtr.vcf.gz:LongTR" \
           "${SAMPLE}.merged.vcf:MERGED VCF"; do
    p="${SAMPLE_OUT}/${f%%:*}"; label="${f#*:}"
    if [[ -f "${p}" ]]; then printf '  [x] %-44s %s\n' "${label}" "$(du -h "${p}" | cut -f1)"
    else printf '  [ ] %-44s\n' "${label}"; fi
  done
  echo
  [[ -f "${LOG}" ]] && { echo "last 5 log lines:"; tail -5 "${LOG}" | sed 's/^/  /'; }
  exit 0
fi

# --- threads ----------------------------------------------------------
if [[ -z "${THREADS}" ]]; then
  NCPU="$(getconf _NPROCESSORS_ONLN 2>/dev/null || nproc 2>/dev/null || echo 4)"
  THREADS=$(( NCPU > 1 ? NCPU - 1 : 1 ))
  [[ "${THREADS}" -gt 16 ]] && THREADS=16
fi

# --- preflight --------------------------------------------------------
preflight() {
  local rc=0
  note "preflight"

  command -v str-toolkit >/dev/null \
    && echo "  [ok] str-toolkit on PATH" \
    || { echo "  [!!] str-toolkit not on PATH. Run: pip install -e ."; rc=1; }

  [[ -f "${CONFIG}" ]] \
    && echo "  [ok] config ${CONFIG}" \
    || { echo "  [!!] config not found: ${CONFIG}"; rc=1; }

  if [[ -n "${BAM}" ]]; then
    [[ -f "${BAM}" ]] && echo "  [ok] bam $(du -h "${BAM}" | cut -f1)" \
      || { echo "  [!!] bam not found: ${BAM}"; rc=1; }
  fi
  if [[ -n "${FASTQ}" ]]; then
    [[ -f "${FASTQ}" ]] && echo "  [ok] fastq $(du -h "${FASTQ}" | cut -f1)" \
      || { echo "  [!!] fastq not found: ${FASTQ}"; rc=1; }
  fi
  [[ -n "${BAM}" || -n "${FASTQ}" ]] || { echo "  [!!] give --bam and/or --fastq"; rc=1; }

  if [[ -n "${REGIONS}" ]]; then
    if [[ -f "${REGIONS}" ]]; then
      echo "  [ok] regions ${REGIONS} ($(wc -l < "${REGIONS}") loci)"
    else
      echo "  [!!] regions file not found: ${REGIONS}"; rc=1
    fi
  fi

  # reference and catalogs named in the config
  if [[ -f "${CONFIG}" ]]; then
    while IFS= read -r line; do
      case "${line}" in
        *:\ /*)
          key="${line%%:*}"; key="${key##*[[:space:]]}"
          val="${line#*: }"; val="${val%%[[:space:]]*}"
          case "${key}" in
            reference|catalog|repeats_bed|regions_bed|mmi)
              if [[ -e "${val}" ]]; then echo "  [ok] ${key} ${val}"
              else echo "  [!!] ${key} missing: ${val}"; rc=1; fi ;;
            model_prefix|last_ref_db)
              # a directory for clair3, a path prefix for lastdb
              if [[ -e "${val}" ]] || compgen -G "${val}*" >/dev/null; then
                echo "  [ok] ${key} ${val}"
              else echo "  [!!] ${key} missing: ${val}"; rc=1; fi ;;
          esac ;;
      esac
    done < "${CONFIG}"
  fi

  # the environment runner, then the binary inside each environment
  RUNNER="$(awk -F': *' '/^env_runner:/ {print $2}' "${CONFIG}" 2>/dev/null | tr -d '[:space:]')"
  if [[ -z "${RUNNER}" ]]; then
    for c in micromamba mamba conda; do
      command -v "$c" >/dev/null && { RUNNER="$c"; break; }
    done
  fi
  if [[ -n "${RUNNER}" ]] && { command -v "${RUNNER}" >/dev/null || [[ -x "${RUNNER}" ]]; }; then
    echo "  [ok] environment runner: ${RUNNER}"
    # The environment and binary names come from config.yaml through the
    # package's own parser, rather than being repeated here where they would
    # drift from the defaults the config actually overrides.
    local spec
    spec="$(python3 - "${CONFIG}" "${TOOLS:-vamos tandem-genotypes longtr}" <<'PYEOF'
import sys
from str_toolkit.config import Config

cfg = Config.from_yaml(sys.argv[1])
tools = sys.argv[2].split()
rows = []
if "vamos" in tools:
    rows += [(cfg.vamos.env_clair3, cfg.vamos.bin_clair3),
             (cfg.vamos.env_whatshap, cfg.vamos.bin_whatshap),
             (cfg.vamos.env_vamos, cfg.vamos.bin_vamos)]
if "longtr" in tools:
    rows += [(cfg.longtr.env, cfg.longtr.bin)]
if "trgt" in tools:
    rows += [(cfg.trgt.env, cfg.trgt.bin)]
if "tandem-genotypes" in tools:
    ctg = cfg.tandem_genotypes
    rows += [(ctg.env_last, ctg.bin_last_train),
             (ctg.env_last, ctg.bin_lastal),
             (ctg.env_last, ctg.bin_last_split),
             (ctg.env_tandem, ctg.bin_tandem_genotypes)]
for env, binary in rows:
    # '|' rather than a tab: bash strips an empty leading field when IFS
    # holds whitespace, which would silently drop every tool configured
    # with no environment (an absolute binary path).
    print(f"{env}|{binary}")
PYEOF
)"
    if [[ -z "${spec}" ]]; then
      echo "  [!!] could not read the tool settings from ${CONFIG}"; rc=1
    else
      while IFS='|' read -r env binary; do
        [[ -n "${binary}" ]] || continue
        check_env "${env}" "${binary}" || rc=1
      done <<<"${spec}"
    fi

    # minimap2 is only needed when an alignment actually runs, which is when
    # no BAM was given or --realign was requested. align_env in the config
    # decides where it is looked for.
    if [[ -z "${BAM}" || "${REALIGN}" -eq 1 ]]; then
      local aenv
      aenv="$(awk -F': *' '/^align_env:/ {print $2}' "${CONFIG}" 2>/dev/null | tr -d '[:space:]')"
      if [[ -n "${aenv}" ]]; then
        check_env "${aenv}" minimap2 || rc=1
        check_env "${aenv}" samtools || rc=1
      elif command -v minimap2 >/dev/null; then
        echo "  [!] minimap2 on PATH ($(command -v minimap2)), no align_env set."
        echo "      It will be looked for inside the genotyping tool's environment"
        echo "      first. Set align_env in ${CONFIG} to avoid a surprise."
      else
        echo "  [!!] an alignment is needed but minimap2 was not found. Set align_env"
        echo "       in ${CONFIG} to an environment that has minimap2 and samtools."
        rc=1
      fi
    else
      echo "  [ok] alignment not needed, the input BAM is reused as is"
    fi
  else
    echo "  [!!] no environment runner. Install micromamba, mamba or conda, or set"
    echo "       env_runner in ${CONFIG} to the full path of one of them."
    rc=1
  fi

  # samtools, needed to index the BAM or to extract reads from it
  local senv
  senv="$(awk -F': *' '/^samtools_env:/ {print $2}' "${CONFIG}" 2>/dev/null)"
  if [[ -n "${senv}" ]]; then
    check_env "${senv}" samtools || rc=1
  elif command -v samtools >/dev/null; then
    echo "  [ok] samtools on PATH"
  else
    echo "  [!!] no samtools. Set samtools_env in ${CONFIG} or put samtools on PATH"; rc=1
  fi

  # free space where the outputs go
  mkdir -p "${OUTDIR}" 2>/dev/null
  local avail
  avail="$(df -BG --output=avail "${OUTDIR}" 2>/dev/null | tail -1 | tr -dc '0-9')"
  if [[ -n "${avail}" ]]; then
    echo "  [ok] free space in ${OUTDIR}: ${avail} G"
    # Only tandem-genotypes needs reads, so the FASTQ extraction and the much
    # larger MAF only cost space when that tool is in the list. LAST writes
    # the MAF as plain text carrying both aligned sequences, so it commonly
    # lands between 1.5 and 3 times the size of the input BAM, and it has to
    # exist in full before tandem-genotypes reads it.
    local want_tg=0
    case " ${TOOLS:-vamos tandem-genotypes longtr} " in *" tandem-genotypes "*) want_tg=1 ;; esac
    if [[ "${want_tg}" -eq 1 && -n "${BAM}" ]]; then
      local bam_g fastq_g maf_lo maf_hi need_lo need_hi
      # stat reports the apparent size, which du understates for a sparse file
      bam_g="$(( $(stat -c %s "${BAM}" 2>/dev/null || echo 0) / 1073741824 ))"
      if [[ -n "${bam_g}" && "${bam_g}" -gt 0 ]]; then
        if [[ -n "${FASTQ}" ]]; then fastq_g=0; else fastq_g=$(( bam_g / 2 )); fi
        maf_lo=$(( bam_g * 3 / 2 )); maf_hi=$(( bam_g * 3 ))
        need_lo=$(( fastq_g + maf_lo )); need_hi=$(( fastq_g + maf_hi ))
        echo "  tandem-genotypes space estimate, from a ${bam_g} G BAM:"
        [[ "${fastq_g}" -gt 0 ]] && echo "    extracted FASTQ      ~${fastq_g} G"
        echo "    LAST MAF             ~${maf_lo} to ${maf_hi} G, needed in full at once"
        echo "    peak requirement     ~${need_lo} to ${need_hi} G against ${avail} G free"
        if [[ "${avail}" -lt "${need_lo}" ]]; then
          echo "  [!!] not enough space for tandem-genotypes on this filesystem."
          echo "       The run would fill the disk partway through lastal, after hours"
          echo "       of work. Either free space, point --outdir at a larger"
          echo "       filesystem, or run the other tools now and tandem-genotypes"
          echo "       on a subsample or on the cluster."
          rc=1
        elif [[ "${avail}" -lt "${need_hi}" ]]; then
          echo "  [!] space is enough at the low end of the estimate only. Watch it"
          echo "      during lastal: df -BG $(dirname "${OUTDIR}")"
        else
          echo "  [ok] space is sufficient across the estimate"
        fi
      fi
    fi
  fi

  echo "  threads: ${THREADS}"
  return "${rc}"
}

check_env() {
  # Reports where a binary actually resolves from inside an environment.
  #
  # A conda environment keeps the system directories on PATH, so plain
  # `command -v` inside it also finds /usr/bin/<bin> and would report the
  # environment as providing a tool it does not contain. The resolved path is
  # therefore compared against the environment prefix, and a hit from outside
  # is reported as such rather than as a pass.
  local env="$1" bin="$2"

  # No environment: the binary is either an absolute path or on the PATH.
  if [[ -z "${env}" ]]; then
    if [[ "${bin}" == /* ]]; then
      if [[ -x "${bin}" ]]; then echo "  [ok] ${bin} (absolute path, no environment)"; return 0
      else echo "  [!!] ${bin} is not an executable file"; return 1; fi
    fi
    if command -v "${bin}" >/dev/null; then
      echo "  [ok] ${bin} on PATH ($(command -v "${bin}")), no environment"
      return 0
    fi
    echo "  [!!] ${bin} not on PATH and no environment configured for it"
    return 1
  fi

  # An absolute path is used as is, so the environment is irrelevant to
  # finding it, though it still supplies the libraries the binary may need.
  if [[ "${bin}" == /* ]]; then
    if [[ -x "${bin}" ]]; then echo "  [ok] ${bin} (absolute path, run inside env ${env})"; return 0
    else echo "  [!!] ${bin} is not an executable file"; return 1; fi
  fi

  local -a pre=("${RUNNER}" run -n "${env}")
  [[ "$(basename "${RUNNER}")" == conda* ]] && pre=("${RUNNER}" run --no-capture-output -n "${env}")

  local resolved prefix
  resolved="$("${pre[@]}" bash -c "command -v ${bin}" 2>/dev/null | tr -d '\r' | tail -1)"
  if [[ -z "${resolved}" ]]; then
    echo "  [!!] ${bin} not found at all from env ${env}"
    return 1
  fi

  prefix="$("${pre[@]}" bash -c 'echo "${CONDA_PREFIX:-}"' 2>/dev/null | tr -d '\r' | tail -1)"
  if [[ -n "${prefix}" && "${resolved}" == "${prefix}"/* ]]; then
    echo "  [ok] env ${env} provides ${bin}"
    return 0
  fi
  # Usable, but it is a system or base install rather than this environment's,
  # so the version is whatever the machine has.
  echo "  [!] ${bin} resolves to ${resolved}, outside env ${env}"
  echo "      It will run, but check its version is the one you intend."
  return 0
}

mkdir -p logs "${OUTDIR}"

if [[ "${PREFLIGHT_ONLY}" -eq 1 ]]; then
  preflight; exit $?
fi

# --- reduced regions catalog -----------------------------------------
# Rewrites only longtr.regions_bed into a copy of the config, so the
# original config.yaml is never touched.
EFFECTIVE_CONFIG="${CONFIG}"
if [[ -n "${REGIONS}" ]]; then
  [[ -f "${REGIONS}" ]] || die "regions file not found: ${REGIONS}"
  EFFECTIVE_CONFIG="${OUTDIR}/.config_${SAMPLE}.yaml"
  sed "s|^\([[:space:]]*\)regions_bed:.*|\1regions_bed: $(readlink -f "${REGIONS}")|" \
    "${CONFIG}" > "${EFFECTIVE_CONFIG}"
  grep -q "$(basename "${REGIONS}")" "${EFFECTIVE_CONFIG}" \
    || die "could not substitute regions_bed into ${CONFIG}, is the key present?"
fi

# --- build the command ------------------------------------------------
CMD=(str-toolkit detect --sample "${SAMPLE}" --config "${EFFECTIVE_CONFIG}"
     --threads "${THREADS}" -o "${OUTDIR}")
[[ -n "${BAM}" ]]   && CMD+=(--bam "${BAM}")
[[ -n "${FASTQ}" ]] && CMD+=(--fastq "${FASTQ}")
[[ -n "${TOOLS}" ]] && { read -r -a _t <<<"${TOOLS}"; CMD+=(--tools "${_t[@]}"); }
[[ "${REALIGN}" -eq 1 ]] && CMD+=(--realign)

if [[ "${DRYRUN}" -eq 1 ]]; then
  printf '%q ' "${CMD[@]}"; echo
  exit 0
fi

# --- detach -----------------------------------------------------------
# Re-exec under setsid + nohup so the run survives the terminal closing.
# _STR_RELAUNCHED marks the second pass, which does the actual work.
if [[ "${BACKGROUND}" -eq 1 && -z "${_STR_RELAUNCHED:-}" ]]; then
  if [[ -f "${PIDFILE}" ]] && kill -0 "$(cat "${PIDFILE}")" 2>/dev/null; then
    die "a run for ${SAMPLE} is already active (pid $(cat "${PIDFILE}")). Use --status."
  fi
  ARGS=(--sample "${SAMPLE}" --config "${CONFIG}" --threads "${THREADS}" --outdir "${OUTDIR}")
  [[ -n "${BAM}" ]]     && ARGS+=(--bam "${BAM}")
  [[ -n "${FASTQ}" ]]   && ARGS+=(--fastq "${FASTQ}")
  [[ -n "${TOOLS}" ]]   && ARGS+=(--tools "${TOOLS}")
  [[ -n "${REGIONS}" ]] && ARGS+=(--regions "${REGIONS}")
  [[ "${REALIGN}" -eq 1 ]] && ARGS+=(--realign)
  _STR_RELAUNCHED=1 setsid nohup "$0" "${ARGS[@]}" >> "${LOG}" 2>&1 < /dev/null &
  echo "started in the background, pid $!"
  echo "log    : ${LOG}      (tail -f ${LOG})"
  echo "status : $0 --sample ${SAMPLE} --outdir ${OUTDIR} --status"
  exit 0
fi

# --- run --------------------------------------------------------------
echo "$$" > "${PIDFILE}"
trap 'rm -f "${PIDFILE}"' EXIT

note "=== ${SAMPLE} on $(hostname), ${THREADS} threads ==="
if ! preflight; then
  note "preflight failed, nothing was run"
  exit 1
fi

note "command: ${CMD[*]}"
# nice and ionice keep a shared machine usable while this runs for hours
RUNNER=(nice -n 10)
command -v ionice >/dev/null && RUNNER=(ionice -c2 -n7 nice -n 10)

START="$(date +%s)"
"${RUNNER[@]}" "${CMD[@]}"
RC=$?
ELAPSED=$(( $(date +%s) - START ))

note "exit code ${RC} after $((ELAPSED / 3600))h$(( (ELAPSED % 3600) / 60 ))m"

if [[ "${RC}" -eq 0 ]]; then
  note "outputs in ${SAMPLE_OUT}"
  ls -lh "${SAMPLE_OUT}" 2>/dev/null | sed 's/^/  /'
  MERGED="${SAMPLE_OUT}/${SAMPLE}.merged.vcf"
  if [[ -f "${MERGED}" ]]; then
    note "merged VCF: $(grep -vc '^#' "${MERGED}") records"
    note "loci by number of supporting tools:"
    grep -v '^#' "${MERGED}" | grep -o 'SOURCES=[^;]*' | sed 's/SOURCES=//' | awk -F',' '
    { delete seen; n = 0
      for (i = 1; i <= NF; i++) { fam = $i; sub(/_[^_]*$/, "", fam)
        if (!(fam in seen)) { seen[fam] = 1; n++ } }
      count[n]++ }
    END { for (k = 1; k <= 4; k++) if (k in count) printf "    %d tool(s): %7d loci\n", k, count[k] }'
  fi
else
  note "the run failed. Every step is idempotent, so relaunching resumes where it"
  note "stopped. Delete ${SAMPLE_OUT} first if you changed the configuration."
fi
exit "${RC}"
