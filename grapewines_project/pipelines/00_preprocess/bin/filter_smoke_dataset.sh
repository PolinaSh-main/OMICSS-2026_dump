#!/bin/bash
#SBATCH --job-name=vcf_smoke
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

set -euo pipefail

WORKDIR="/path/to/workdir/grapewines_project"
SMOKE_DIR="${WORKDIR}/smoke_test"
LOG="${SMOKE_DIR}/logs/smoke_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "${SMOKE_DIR}/logs"

exec > >(tee -a "${LOG}") 2>&1

echo "=== Smoke filtering started: $(date) ==="
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Input:  ${SMOKE_DIR}/smoke_input.vcf.gz"
echo ""

cd "${SMOKE_DIR}"

vcftools \
  --gzvcf "${SMOKE_DIR}/smoke_input.vcf.gz" \
  --remove-filtered-all \
  --min-alleles 2 --max-alleles 2 \
  --max-missing 0.6 \
  --maf 0.005 \
  --recode --recode-INFO-all \
  --out cauc_smoke_filtered

echo ""
echo "Variants after filtering: $(grep -c '^[^#]' cauc_smoke_filtered.recode.vcf)"
echo "Output: ${SMOKE_DIR}/cauc_smoke_filtered.recode.vcf"
echo "=== Smoke filtering finished: $(date) ==="