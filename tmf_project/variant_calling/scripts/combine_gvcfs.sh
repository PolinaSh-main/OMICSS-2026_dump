#!/bin/bash
set -euo pipefail

REF=$1
OUTPUT=$2
GATK=$3

shift 3

GVCFS=("$@")


for gvcf in "${GVCFS[@]}"
do
    if [ ! -f "$gvcf" ]; then
        echo "Missing $gvcf"
        exit 1
    fi
done


$GATK --java-options "-Xmx6g" CombineGVCFs \
    -R "$REF" \
    $(printf -- "--variant %s " "${GVCFS[@]}") \
    -O "$OUTPUT"


test -f "$OUTPUT"