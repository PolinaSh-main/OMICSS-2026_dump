#!/bin/bash
set -euo pipefail

REF=$1
INPUT=$2
OUTPUT=$3
GATK=$4


test -f "$INPUT"


$GATK --java-options "-Xmx6g" GenotypeGVCFs \
    -R "$REF" \
    -V "$INPUT" \
    -O "$OUTPUT"


test -f "$OUTPUT"