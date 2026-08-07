#!/bin/bash
set -euo pipefail

INPUT=$1
OUTPUT=$2
GATK=$3

if [ ! -f "$INPUT" ]; then
    echo "ERROR: missing input BAM: $INPUT"
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT")"

$GATK SortSam \
    -I "$INPUT" \
    -O "$OUTPUT" \
    --SORT_ORDER coordinate

test -f "$OUTPUT"

echo "SortSam finished"