#!/bin/bash
#SBATCH --mem=10gb
#SBATCH --cpus-per-task=4
#SBATCH --partition=thin
#SBATCH --output=log/slurm_fastp-%x-%j.log
#SBATCH -o log/admixturegp3.out
#SBATCH -e log/admixturegp3.err

#to run:
#sbatch -D /mnt/nas1/proj/admixture/results/polina_sh/logomicss26/gp3/admixture/results/polina_sh \
#    /mnt/nas1/proj/omicss26/gp3/admixture/script/admixture_polia.sh

# Please doubple check the paths before running the script!
admixture="/mnt/nas1/proj/omicss26/soft/admixture_linux-1.3.0/admixture"
data="/mnt/nas1/proj/omicss26/gp3/admixture/data_filtered/cauc_filtered.final.bed"
log_dir="/mnt/nas1/proj/omicss26/gp3/admixture/results/polina_sh/log"

echo "creating lod directory"
mkdir -p "$log_dir"

echo "Running the admixture ..."

for K in {2..10}; do
    echo "=== Running ADMIXTURE for K=$K ==="
    ${admixture} --cv ${data} $K -j4 | tee "${log_dir}/log${K}.out"
done

echo "Admixture done!"