#!/bin/bash
set -euo pipefail

INPUT=$1
OUTPUT=$2
METRICS=$3
GATK=$4

if [ ! -f "$INPUT" ]; then
    echo "ERROR: missing input BAM"
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT")"

$GATK --java-options "-Xmx6g" MarkDuplicates \
    -I "$INPUT" \
    -O "$OUTPUT" \
    -M "$METRICS" \
    --CREATE_INDEX true


test -f "$OUTPUT"

echo "MarkDuplicates finished"