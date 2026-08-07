#!/bin/bash
#SBATCH --job-name=enrichment
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=log/enrichment_%j.out
#SBATCH --error=log/enrichment_%j.err

python enrichment_analysis.py