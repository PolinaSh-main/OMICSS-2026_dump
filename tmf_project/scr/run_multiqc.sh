#!/bin/bash

set -euo pipefail

PROJECT_DIR=$(dirname "$(dirname "$(realpath "$0")")")

FASTQC_DIR="$PROJECT_DIR/results/fastqc"
MULTIQC_DIR="$PROJECT_DIR/results/multiqc"
LOG_DIR="$PROJECT_DIR/logs"

mkdir -p "$MULTIQC_DIR"
mkdir -p "$LOG_DIR"


echo "Running MultiQC..."

multiqc \
    "$FASTQC_DIR" \
    --outdir "$MULTIQC_DIR" \
    > "$LOG_DIR/multiqc.log" 2>&1


echo "MultiQC finished successfully"
echo "Report: $MULTIQC_DIR/multiqc_report.html"