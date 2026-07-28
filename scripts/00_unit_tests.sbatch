#!/bin/bash
#SBATCH --job-name=str-unit-tests
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --partition=<ta_partition>

# Ne teste QUE la logique Python (parsing, clustering, CLI) -- ne lance
# aucun des 3 outils bioinfo, donc pas besoin d'environnements micromamba ici.
# À lancer en tout premier pour valider que le package est bien installé.

set -euo pipefail
mkdir -p logs

cd "${SLURM_SUBMIT_DIR}"
python -m pytest tests/ -v
