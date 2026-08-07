#!/bin/bash
#SBATCH --job-name=bwa_index
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/bwa_index_%j.out
#SBATCH --error=logs/bwa_index_%j.err

mkdir -p logs

cd /mnt/nas0/user/polina.shevyakova/tmf_project/data/hg38_index

bwa index hg38.fa