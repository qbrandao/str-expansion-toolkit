#!/bin/bash
# Locate the tool binaries across the conda-style environments on this
# machine, and print the config.yaml block that matches what was found.
#
# The environment names in config.example.yaml (clair3, whatshap-env, vamos,
# longtr, last_env, tandem-env) are the ones used on the machine this project
# was developed on. Another machine names them differently, or installs a
# tool outside any environment, and then the preflight reports the binary as
# missing even though it is present. This script says where each one actually
# is.
#
# Usage:
#   ./scripts/find_tools.sh [runner]
#
# runner is micromamba, mamba, conda, or a full path to one of them. Omitted,
# the first one found on PATH is used.

set -uo pipefail

RUNNER="${1:-}"
if [[ -z "${RUNNER}" ]]; then
  for c in micromamba mamba conda; do
    command -v "$c" >/dev/null && { RUNNER="$c"; break; }
  done
fi
[[ -n "${RUNNER}" ]] || { echo "ERROR: no micromamba, mamba or conda found" >&2; exit 2; }
echo "runner: ${RUNNER}"
echo

# --- environment inventory -------------------------------------------
# `env list` prints a comment header then one line per environment, whose
# last field is the prefix. The active one carries a '*', which is why the
# path is taken from the last field rather than from a fixed column.
ENV_PATHS=()
while IFS= read -r line; do
  ENV_PATHS+=("${line}")
done < <("${RUNNER}" env list 2>/dev/null | awk '!/^#/ && NF {print $NF}' | grep '^/' | sort -u)

if [[ "${#ENV_PATHS[@]}" -eq 0 ]]; then
  echo "WARNING: '${RUNNER} env list' returned no environment. Only PATH will be searched." >&2
fi

echo "environments found: ${#ENV_PATHS[@]}"
for p in "${ENV_PATHS[@]}"; do
  name="$(basename "${p}")"
  # the base environment's prefix has no /envs/ component
  [[ "${p}" == */envs/* ]] || name="base"
  printf '  %-22s %s\n' "${name}" "${p}"
done
echo

# --- binaries the pipeline needs --------------------------------------
# One line per binary: "binary|config key it belongs to"
REQUIRED=(
  "run_clair3.sh|vamos.env_clair3"
  "whatshap|vamos.env_whatshap"
  "vamos|vamos.env_vamos"
  "samtools|samtools_env"
  "minimap2|longtr.env (alignment)"
  "LongTR|longtr.env"
  "lastal|tandem_genotypes.env_last"
  "last-train|tandem_genotypes.env_last"
  "last-split|tandem_genotypes.env_last"
  "tandem-genotypes|tandem_genotypes.env_tandem"
  "trgt|trgt.env (opt-in)"
)

declare -A FOUND_IN   # binary -> space-separated env names
MISSING=()

echo "----------------------------------------------------------------"
printf '%-18s %-34s %s\n' "BINARY" "CONFIG KEY" "FOUND IN"
echo "----------------------------------------------------------------"

for entry in "${REQUIRED[@]}"; do
  bin="${entry%%|*}"; key="${entry#*|}"
  hits=()
  for p in "${ENV_PATHS[@]}"; do
    if [[ -x "${p}/bin/${bin}" ]]; then
      name="$(basename "${p}")"
      # the base environment's prefix has no /envs/ component
      [[ "${p}" == */envs/* ]] || name="base"
      hits+=("${name}")
    fi
  done
  # also the current PATH, for a tool installed outside any environment
  if command -v "${bin}" >/dev/null 2>&1; then
    hits+=("PATH:$(command -v "${bin}")")
  fi

  if [[ "${#hits[@]}" -eq 0 ]]; then
    printf '%-18s %-34s %s\n' "${bin}" "${key}" "NOT FOUND"
    MISSING+=("${bin}")
  else
    printf '%-18s %-34s %s\n' "${bin}" "${key}" "${hits[*]}"
    FOUND_IN["${bin}"]="${hits[*]}"
  fi
done
echo "----------------------------------------------------------------"
echo

# --- suggested config -------------------------------------------------
# A binary often sits in several environments, so the pick is not simply the
# first hit: an environment named after the tool is almost always the one
# installed for it (whatshap in 'whatshap' rather than in 'clair3', which
# merely has it as a dependency). Falling back to the first environment hit,
# and to an empty answer when the only hit is the PATH, which the config
# expresses by leaving the key out.
pick() {
  local bin="$1" hits lowbin
  hits="${FOUND_IN[${bin}]:-}"
  [[ -n "${hits}" ]] || { echo ""; return; }
  lowbin="$(tr '[:upper:]' '[:lower:]' <<<"${bin}")"

  # an environment whose name contains the binary name
  for h in ${hits}; do
    case "${h}" in PATH:*) continue ;; esac
    case "$(tr '[:upper:]' '[:lower:]' <<<"${h}")" in
      *"${lowbin}"*) echo "${h}"; return ;;
    esac
  done
  # otherwise the first real environment
  for h in ${hits}; do
    case "${h}" in PATH:*) continue ;; esac
    echo "${h}"; return
  done
  echo ""
}

on_path() {
  case " ${FOUND_IN[$1]:-} " in *" PATH:"*) return 0 ;; esac
  return 1
}

E_CLAIR3="$(pick run_clair3.sh)"
E_WHATSHAP="$(pick whatshap)"
E_VAMOS="$(pick vamos)"
E_SAMTOOLS="$(pick samtools)"
E_LONGTR="$(pick LongTR)"
E_LAST="$(pick lastal)"
E_TANDEM="$(pick tandem-genotypes)"
E_TRGT="$(pick trgt)"

echo "Suggested config.yaml environment settings, based on what was found:"
echo
echo "  env_runner: $(command -v "${RUNNER}" 2>/dev/null || echo "${RUNNER}")"

# samtools and minimap2 are usually system installs. When one is on the PATH
# the key is better left out, which keeps the version explicit and avoids
# pinning the pipeline to an unrelated environment.
if on_path samtools; then
  echo "  # samtools_env is left out on purpose: samtools is on the PATH at"
  echo "  #   $(command -v samtools). Check 'samtools --version' is recent enough."
elif [[ -n "${E_SAMTOOLS}" ]]; then
  echo "  samtools_env: ${E_SAMTOOLS}"
else
  echo "  # samtools_env: <samtools NOT FOUND>"
fi

E_MINIMAP="$(pick minimap2)"
if on_path minimap2; then
  echo "  # align_env is left out on purpose: minimap2 is on the PATH at"
  echo "  #   $(command -v minimap2). It only matters with --realign or a FASTQ-only run."
elif [[ -n "${E_MINIMAP}" ]]; then
  echo "  align_env: ${E_MINIMAP}   # needs BOTH minimap2 and samtools"
else
  echo "  # align_env: <minimap2 NOT FOUND, only needed with --realign or a FASTQ-only run>"
fi
echo
echo "  vamos:"
echo "    env_clair3: ${E_CLAIR3:-<NOT FOUND>}"
echo "    env_whatshap: ${E_WHATSHAP:-<NOT FOUND>}"
echo "    env_vamos: ${E_VAMOS:-<NOT FOUND>}"
echo
echo "  longtr:"
echo "    env: ${E_LONGTR:-<NOT FOUND>}"
echo
echo "  tandem_genotypes:"
echo "    env_last: ${E_LAST:-<NOT FOUND>}"
echo "    env_tandem: ${E_TANDEM:-<NOT FOUND>}"
echo
echo "  trgt:"
echo "    env: ${E_TRGT:-<NOT FOUND>}   # opt-in, only needed with --tools ... trgt"
echo

if [[ "${#MISSING[@]}" -gt 0 ]]; then
  echo "Not found in any environment or on the PATH: ${MISSING[*]}"
  echo
  echo "Looking for a hand-built copy under \$HOME, /opt and /usr/local:"
  declare -A HANDBUILT=()
  for m in "${MISSING[@]}"; do
    hit="$(find "${HOME}" /opt /usr/local -maxdepth 6 -type f -name "${m}" -perm -u+x 2>/dev/null | head -1)"
    if [[ -n "${hit}" ]]; then
      HANDBUILT["${m}"]="${hit}"
      printf '  %-18s %s\n' "${m}" "${hit}"
    else
      printf '  %-18s not found\n' "${m}"
    fi
  done
  echo

  if [[ "${#HANDBUILT[@]}" -gt 0 ]]; then
    echo "Those were built outside any environment. The config expresses that with"
    echo "an absolute bin path and an empty env name, so add this:"
    echo
    for m in "${!HANDBUILT[@]}"; do
      case "${m}" in
        vamos)            echo "  vamos:"
                          echo "    env_vamos: \"\""
                          echo "    bin_vamos: ${HANDBUILT[$m]}" ;;
        run_clair3.sh)    echo "  vamos:"
                          echo "    env_clair3: \"\""
                          echo "    bin_clair3: ${HANDBUILT[$m]}" ;;
        whatshap)         echo "  vamos:"
                          echo "    env_whatshap: \"\""
                          echo "    bin_whatshap: ${HANDBUILT[$m]}" ;;
        LongTR)           echo "  longtr:"
                          echo "    env: \"\""
                          echo "    bin: ${HANDBUILT[$m]}" ;;
        trgt)             echo "  trgt:"
                          echo "    env: \"\""
                          echo "    bin: ${HANDBUILT[$m]}" ;;
        lastal)           echo "  tandem_genotypes:"
                          echo "    bin_lastal: ${HANDBUILT[$m]}" ;;
        last-train)       echo "  tandem_genotypes:"
                          echo "    bin_last_train: ${HANDBUILT[$m]}" ;;
        last-split)       echo "  tandem_genotypes:"
                          echo "    bin_last_split: ${HANDBUILT[$m]}" ;;
        tandem-genotypes) echo "  tandem_genotypes:"
                          echo "    bin_tandem_genotypes: ${HANDBUILT[$m]}" ;;
        minimap2|samtools)
                          echo "  # ${m} at ${HANDBUILT[$m]} is not reached through a config key."
                          echo "  # Put its directory on PATH, or set align_env / samtools_env to an"
                          echo "  # environment that has it." ;;
      esac
    done
    echo
    echo "Merge the keys under one block per section rather than repeating the"
    echo "section header. Keeping the env name set alongside an absolute path also"
    echo "works, and is worth doing when the binary needs the shared libraries of"
    echo "that environment."
  fi

  still_missing=()
  for m in "${MISSING[@]}"; do
    [[ -n "${HANDBUILT[$m]:-}" ]] || still_missing+=("${m}")
  done
  if [[ "${#still_missing[@]}" -gt 0 ]]; then
    echo
    echo "Genuinely absent from this machine: ${still_missing[*]}"
    echo "Install into an environment, or build and point the matching bin_* key"
    echo "at the result."
  fi
  echo
  echo "Which tools you actually need depends on what you run:"
  echo "  LongTR alone            LongTR, minimap2, samtools"
  echo "  tandem-genotypes alone  lastal, last-train, last-split, tandem-genotypes, samtools"
  echo "  VAMOS alone             run_clair3.sh, whatshap, vamos, samtools"
  echo "So a missing lastal does not stop a LongTR run."
fi
