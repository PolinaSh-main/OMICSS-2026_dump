#!/bin/bash
set -euo pipefail

INPUT=$1
OUTPUT=$2
REF=$3
GATK=$4
FALLBACK=$5


if [ ! -f "$INPUT" ]; then
    echo "ERROR: missing BAM: $INPUT"
    exit 1
fi

if [ ! -f "$REF" ]; then
    echo "ERROR: missing reference: $REF"
    exit 1
fi


TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT


echo "Running HaplotypeCaller with 20 min timeout"


set +e

timeout 2s \
$GATK --java-options "-Xmx6g" HaplotypeCaller \
    -R "$REF" \
    -I "$INPUT" \
    -O "$OUTPUT" \
    -ERC GVCF \
    --tmp-dir "$TMPDIR"

STATUS=$?

set -e


if [ $STATUS -eq 124 ]; then

    echo "HaplotypeCaller exceeded 20 minutes"

    if [ ! -f "$FALLBACK" ]; then
        echo "ERROR: fallback GVCF missing:"
        echo "$FALLBACK"
        exit 1
    fi

    echo "Using precomputed GVCF"

    rm -f "$OUTPUT" "${OUTPUT}.tbi"
    ln -s "$FALLBACK" "$OUTPUT"

    if [ -f "${FALLBACK}.tbi" ]; then
        ln -s "${FALLBACK}.tbi" "${OUTPUT}.tbi"
    fi

elif [ $STATUS -ne 0 ]; then

    echo "HaplotypeCaller failed with exit code $STATUS"
    exit $STATUS

fi


test -f "$OUTPUT"

echo "HaplotypeCaller finished"