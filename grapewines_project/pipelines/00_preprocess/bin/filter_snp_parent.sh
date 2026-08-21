#!/bin/bash
#SBATCH --mem=10gb
#SBATCH --partition=thin
#SBATCH --cpus-per-task=4            
#SBATCH --job-name=preprocessing
#SBATCH --output=/path/to/shared-data/admixture/log/filtering_00%j.log  # %j will be replaced with the job ID

input_file=/path/to/shared-data/cauc_filtered.recode.vcf
store_plink=/path/to/shared-data/data/plink
output_dir=/path/to/shared-data/data/plink/final_filtered

script/filter_snp_child.sh $input_file $store_plink $ids $output_dir
