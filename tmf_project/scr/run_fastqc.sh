#!/bin/bash

set -euo pipefail

# Корень проекта: на два уровня выше src/run_fastqc.sh
PROJECT_DIR=$(dirname "$(dirname "$(realpath "$0")")")

# Пути
DATA_DIR="$PROJECT_DIR/data/fastq"
RESULTS_DIR="$PROJECT_DIR/results/fastqc"
LOG_DIR="$PROJECT_DIR/logs"


# Создаем папки, если их нет
mkdir -p "$RESULTS_DIR"
mkdir -p "$LOG_DIR"


echo "Project directory: $PROJECT_DIR"
echo "Input data: $DATA_DIR"
echo "Output: $RESULTS_DIR"


echo "Running FastQC..."

fastqc \
    "$DATA_DIR"/*.fq.gz \
    --outdir "$RESULTS_DIR" \
    --threads 4 \
    > "$LOG_DIR/fastqc.log" 2>&1


echo "FastQC finished successfully"