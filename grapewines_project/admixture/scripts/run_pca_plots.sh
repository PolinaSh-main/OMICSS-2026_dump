#!/bin/bash
#SBATCH --mem=10gb
#SBATCH --cpus-per-task=1
#SBATCH --partition=thin
#SBATCH -o /mnt/nas0/user/polina.shevyakova/grapewines_project/admixture/log/pca.out
#SBATCH -e /mnt/nas0/user/polina.shevyakova/grapewines_project/admixture/log/pca.err

PROJECT="/mnt/nas0/user/polina.shevyakova/grapewines_project/admixture"

mkdir -p "$PROJECT/log"
mkdir -p "$PROJECT/pca_results"

plink \
    --bfile "$PROJECT/data_filtered/cauc_filtered.final" \
    --pca 20 \
    --out "$PROJECT/pca_results/cauc_pca"

python3 "/mnt/nas0/user/polina.shevyakova/grapewines_project/admixture/scripts/pca_plot.py" \
    --eigenvec "$PROJECT/pca_results/cauc_pca.eigenvec" \
    --eigenval "$PROJECT/pca_results/cauc_pca.eigenval" \
    --metadata "/mnt/nas1/proj/omicss26/gp3/data/metadata/cauc_grape_metadata.csv" \
    --outdir "$PROJECT/pca_results"