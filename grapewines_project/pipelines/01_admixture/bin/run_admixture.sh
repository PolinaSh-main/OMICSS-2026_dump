#!/bin/bash
#SBATCH --mem=10gb
#SBATCH --cpus-per-task=2
#SBATCH --partition=thin
#SBATCH --output=log/slurm_fastp-%x-%j.log
#SBATCH -o log/admixturegp3.out
#SBATCH -e log/admixturegp3.err

# Please double check the paths before running the script!
admixture=""/mnt/nas0/proj/vine/shared_files/soft/admixture_linux-1.3.0/admixture""
data="../data/plink/final_filtered/cauc_filtered.final.bed"
log_dir="log"

echo "creating lod directory"
mkdir -p "$log_dir"

echo "Running the admixture ..."

for K in {2..10}; do
    echo "=== Running ADMIXTURE for K=$K ==="
    ${admixture} --cv ${data} $K -j4 | tee "${log_dir}/log${K}.out"
done

echo "Admixture done!"