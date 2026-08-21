#!/bin/bash
#SBATCH --job-name=gatk_sort
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/gatk_sort_%j.out
#SBATCH --error=logs/gatk_sort_%j.err

sample=$1

data_dir='/path/to/workdir/tmf_project/variant_calling'

ref='/path/to/shared-data/ngs_data_analysis/alignment_samtools/ref_genome/hg38.fa'

gatk_bin="/path/to/tools/gatk-4.2.6.1/gatk"

${gatk_bin} SortSam \
  -I ${data_dir}/bam/${sample}_sorted.bam \
  -O ${data_dir}/bam_clean/${sample}_sorted.bam \
  --SORT_ORDER coordinate

${gatk_bin} MarkDuplicates \
  -I ${data_dir}/bam_clean/${sample}_sorted.bam \
  -O ${data_dir}/bam_clean/${sample}_dedup.bam \
  -M ${data_dir}/bam_clean/${sample}_dedup_metrics.txt \
  --CREATE_INDEX true

gvcf_out=${data_dir}/gvcf/${sample}.g.vcf.gz

timeout 20m ${gatk_bin} --java-options "-Xmx6g" HaplotypeCaller \
  -R ${ref} \
  -I ${data_dir}/bam_clean/${sample}_dedup.bam \
  -O ${gvcf_out} \
  -ERC GVCF

if [ $? -eq 124 ]; then
    echo "HaplotypeCaller exceeded 20 minutes. Copying precomputed GVCF..."

    cp /path/to/shared-data/ngs_data_analysis/variant_calling/data/gvcf/${sample}.g.vcf.gz \
       ${gvcf_out}

    cp /path/to/shared-data/ngs_data_analysis/variant_calling/data/gvcf/${sample}.g.vcf.gz.tbi \
       ${gvcf_out}.tbi 2>/dev/null || true
fi

${gatk_bin} CombineGVCFs \
  -R ${ref} \
  --variant ${data_dir}/gvcf/wes46.g.vcf.gz \
  --variant ${data_dir}/gvcf/wes78.g.vcf.gz \
  -O ${data_dir}/gvcf/combined.g.vcf.gz

${gatk_bin} GenotypeGVCFs \
  -R ${ref} \
  -V ${data_dir}/gvcf/combined.g.vcf.gz \
  -O ${data_dir}/vcf/genotyped_variants.vcf.gz

${gatk_bin} SelectVariants \
  -R ${ref} \
  -V ${data_dir}/vcf/genotyped_variants.vcf.gz \
  --select-type-to-include SNP \
  -O ${data_dir}/vcf/snp_variants.vcf

${gatk_bin} SelectVariants \
  -R ${ref} \
  -V ${data_dir}/vcf/genotyped_variants.vcf.gz \
  --select-type-to-include INDEL \
  -O ${data_dir}/vcf/indel_variants.vcf

${gatk_bin} VariantFiltration \
  -R ${ref} \
  -V ${data_dir}/vcf/snp_variants.vcf \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 10.0" --filter-name "SOR10" \
  -filter "FS > 60.0" --filter-name "FS60" \
  -O ${data_dir}/vcf/filtered_snps.vcf

${gatk_bin} VariantFiltration \
  -R ${ref} \
  -V ${data_dir}/vcf/indel_variants.vcf \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 10.0" --filter-name "SOR10" \
  -filter "FS > 100.0" --filter-name "FS100" \
  -O ${data_dir}/vcf/filtered_indels.vcf
