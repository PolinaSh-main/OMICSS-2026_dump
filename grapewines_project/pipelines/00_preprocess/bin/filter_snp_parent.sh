#!/bin/bash
#SBATCH --mem=10gb
#SBATCH --partition=thin
#SBATCH --cpus-per-task=4            
#SBATCH --job-name=preprocessing
#SBATCH --output=/mnt/nas1/proj/omicss26/gp3/admixture/log/filtering_00%j.log  # %j will be replaced with the job ID

input_file=/mnt/nas1/proj/omicss26/gp3/cauc_filtered.recode.vcf
store_plink=/mnt/nas1/proj/omicss26/gp3/data/plink
output_dir=/mnt/nas1/proj/omicss26/gp3/data/plink/final_filtered

script/filter_snp_child.sh $input_file $store_plink $ids $output_dir
