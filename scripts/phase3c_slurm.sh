#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=phase3c.%j.out
set -euo pipefail
python scripts/build_pathogen_db.py --manifest "$1" --genomes-dir "$2" --build-skani --threads "${SLURM_CPUS_PER_TASK:-16}"
