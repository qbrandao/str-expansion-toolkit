#!/bin/bash
#SBATCH --job-name=str-unit-tests
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --partition=<your_partition>

# Tests ONLY the Python logic (parsing, clustering, CLI) -- does not run
# any of the bioinformatics tools, so no micromamba environments are
# needed here. Run this first to validate the package is properly installed.

set -euo pipefail
mkdir -p logs

cd "${SLURM_SUBMIT_DIR}"
python -m pytest tests/ -v
