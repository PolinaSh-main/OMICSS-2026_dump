#!/bin/bash
#SBATCH --mem 10gb
#SBATCH --cpus-per-task=1
#SBATCH --partition=thin
#SBATCH --output=log/slurm_fastp-%x-%j.log
#SBATCH -o log/plinkg_launch.out
#SBATCH -e log/plinkg_launch.err

# Please double check the directories before running the script !

data_dir="../data/vcf"

echo "Converting vcf to plink ..."

plink --vcf ${data_dir}/cauc_filtered.final.vcf.gz --make-bed \
--double-id \
--allow-extra-chr 0 \
--chr-set 19 \
--out ${data_dir}/cauc_filtered.final

echo "Done!"