#!/bin/bash
#SBATCH --job-name=vcf_filter
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G

set -euo pipefail

WORKDIR="/path/to/workdir/grapewines_project"
INPUT_VCF="/path/to/shared-data/data/vcf/cauca_grape.subset.vcf.gz"
LOG="${WORKDIR}/logs/filter_full_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "${WORKDIR}/logs"

exec > >(tee -a "${LOG}") 2>&1

echo "=== Full filtering started: $(date) ==="
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Input:  ${INPUT_VCF}"
echo ""

cd "${WORKDIR}"

vcftools \
  --gzvcf "${INPUT_VCF}" \
  --remove-filtered-all \
  --min-alleles 2 --max-alleles 2 \
  --max-missing 0.6 \
  --maf 0.005 \
  --recode --recode-INFO-all \
  --out cauc_filtered

echo ""
echo "Filtering done. Compressing and indexing ..."

bgzip cauc_filtered.recode.vcf
tabix -p vcf cauc_filtered.recode.vcf.gz

echo "Variants after filtering: $(bcftools view -H cauc_filtered.recode.vcf.gz | wc -l)"
echo "Output: ${WORKDIR}/cauc_filtered.recode.vcf.gz"
echo "=== Full filtering finished: $(date) ==="