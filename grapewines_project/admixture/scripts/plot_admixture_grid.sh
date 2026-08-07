#!/bin/bash
#SBATCH --mem=10gb
#SBATCH --cpus-per-task=1
#SBATCH --partition=thin
#SBATCH -o /mnt/nas0/user/polina.shevyakova/grapewines_project/admixture/log/adm.out
#SBATCH -e /mnt/nas0/user/polina.shevyakova/grapewines_project/admixture/log/adm.err

python plot_admixture_grid.py \
  --qdir   /mnt/nas0/user/polina.shevyakova/grapewines_project/admixture/results/ \
  --orders /mnt/nas0/user/polina.shevyakova/grapewines_project/admixture/output/orders/ \
  --fam    /mnt/nas0/user/polina.shevyakova/grapewines_project/admixture/data_filtered/cauc_filtered.final.fam \
  --metadata /mnt/nas1/proj/omicss26/gp3/data/metadata/cauc_grape_metadata.csv \
  --outdir /mnt/nas0/user/polina.shevyakova/grapewines_project/admixture/output/plots/ \
  --cv     /mnt/nas0/user/polina.shevyakova/grapewines_project/admixture/output/cv/cv_errors.tsv \
  --best-k 7