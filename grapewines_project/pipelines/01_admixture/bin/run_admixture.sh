#!/bin/bash
#SBATCH --mem=10gb
#SBATCH --cpus-per-task=4
#SBATCH --partition=thin
#SBATCH --output=log/slurm_fastp-%x-%j.log
#SBATCH -o log/admixture.out
#SBATCH -e log/admixture.err

#to run:
#sbatch -D /path/to/admixture-run \
#    /path/to/admixture-run/run_admixture.sh

# Please double check the paths before running the script!
admixture="/path/to/tools/admixture_linux-1.3.0/admixture"
data="/path/to/shared-data/admixture/data_filtered/cauc_filtered.final.bed"
log_dir="/path/to/admixture-run/log"

echo "creating lod directory"
mkdir -p "$log_dir"

echo "Running the admixture ..."

for K in {2..10}; do
    echo "=== Running ADMIXTURE for K=$K ==="
    ${admixture} --cv ${data} $K -j4 | tee "${log_dir}/log${K}.out"
done

echo "Admixture done!"