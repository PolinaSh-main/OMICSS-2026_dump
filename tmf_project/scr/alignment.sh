#!/bin/bash
#SBATCH --job-name=alignment
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb
#SBATCH --output=logs/alignment_%j.out
#SBATCH --error=logs/alignment_%j.err


# the sample name is taken from the command line, e.g. sbatch alignment.sh wes46
sample=$1

data_dir='/mnt/nas0/user/polina.shevyakova/tmf_project/data/fastq'
ref='/mnt/nas1/proj/omicss26/ngs_data_analysis/alignment_samtools/ref_genome/hg38.fa'

READ1="${data_dir}/${sample}_chr21_chr16_R1.fastq"
READ2="${data_dir}/${sample}_chr21_chr16_R2.fastq"

mkdir -p bam

bwa mem -t "${SLURM_CPUS_PER_TASK}" \
  -R "@RG\tID:${sample}\tLB:${sample}\tSM:${sample}\tPL:ILLUMINA" \
  "${ref}" "${READ1}" "${READ2}" | \
  samtools view -b - > "bam/${sample}.bam"