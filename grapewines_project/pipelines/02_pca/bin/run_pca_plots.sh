#!/bin/bash
#SBATCH --mem=10gb
#SBATCH --cpus-per-task=1
#SBATCH --partition=thin
#SBATCH -o /path/to/workdir/grapewines_project/admixture/log/pca.out
#SBATCH -e /path/to/workdir/grapewines_project/admixture/log/pca.err

PROJECT="/path/to/workdir/grapewines_project/admixture"

mkdir -p "$PROJECT/log"
mkdir -p "$PROJECT/pca_results"

plink \
    --bfile "$PROJECT/data_filtered/cauc_filtered.final" \
    --pca 20 \
    --out "$PROJECT/pca_results/cauc_pca"

python3 "/path/to/workdir/grapewines_project/admixture/scripts/pca_plot.py" \
    --eigenvec "$PROJECT/pca_results/cauc_pca.eigenvec" \
    --eigenval "$PROJECT/pca_results/cauc_pca.eigenval" \
    --metadata "/path/to/shared-data/data/metadata/cauc_grape_metadata.csv" \
    --outdir "$PROJECT/pca_results"