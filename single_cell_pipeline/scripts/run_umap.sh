#!/bin/bash
#SBATCH --job-name=umap
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=../logs/umap_%j.out
#SBATCH --error=../logs/umap_%j.err

set -euo pipefail

source ../.venv/bin/activate

python run_umap.py \
    --input ../results/integration/harmony/harmony.h5ad \
    --output ../results/umap