#!/bin/bash
set -euo pipefail

REF=$1
INPUT=$2
OUTDIR=$3
GATK=$4


mkdir -p "$OUTDIR"


$GATK SelectVariants \
    -R "$REF" \
    -V "$INPUT" \
    --select-type-to-include SNP \
    -O "$OUTDIR/snps.vcf"


$GATK SelectVariants \
    -R "$REF" \
    -V "$INPUT" \
    --select-type-to-include INDEL \
    -O "$OUTDIR/indels.vcf"



$GATK VariantFiltration \
    -R "$REF" \
    -V "$OUTDIR/snps.vcf" \
    --filter "QD < 2.0" \
    --filter-name "QD2" \
    --filter "QUAL < 30.0" \
    --filter-name "QUAL30" \
    --filter "SOR > 10.0" \
    --filter-name "SOR10" \
    --filter "FS > 60.0" \
    --filter-name "FS60" \
    -O "$OUTDIR/filtered_snps.vcf"



$GATK VariantFiltration \
    -R "$REF" \
    -V "$OUTDIR/indels.vcf" \
    --filter "QD < 2.0" \
    --filter-name "QD2" \
    --filter "QUAL < 30.0" \
    --filter-name "QUAL30" \
    --filter "SOR > 10.0" \
    --filter-name "SOR10" \
    --filter "FS > 100.0" \
    --filter-name "FS100" \
    -O "$OUTDIR/filtered_indels.vcf"