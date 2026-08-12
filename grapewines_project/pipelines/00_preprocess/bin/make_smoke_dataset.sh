#!/bin/bash
# Создаёт маленький subset VCF (первые 2 Мб хромосомы 1) для smoke test

set -euo pipefail

WORKDIR="/mnt/nas0/user/polina.shevyakova/grapewines_project"
INPUT_VCF="/mnt/nas1/proj/omicss26/gp3/data/vcf/cauca_grape.subset.vcf.gz"
SMOKE_DIR="${WORKDIR}/smoke_test"

mkdir -p "${SMOKE_DIR}"

echo "Extracting region 1:1-2000000 ..."
bcftools view \
  -r 1:1-2000000 \
  "${INPUT_VCF}" \
  -Oz -o "${SMOKE_DIR}/smoke_input.vcf.gz"

echo "Indexing ..."
tabix -p vcf "${SMOKE_DIR}/smoke_input.vcf.gz"

N_VARIANTS=$(bcftools view -H "${SMOKE_DIR}/smoke_input.vcf.gz" | wc -l)
echo "Done. Variants in subset: ${N_VARIANTS}"
echo "Output: ${SMOKE_DIR}/smoke_input.vcf.gz"