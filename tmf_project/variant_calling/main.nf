nextflow.enable.dsl = 2

def logDir = "${params.results}/${workflow.runName}/logs"


process SORT_BAM {
    tag "$sample"

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}_sorted.bam")

    script:
    """
    mkdir -p ${logDir}

    (
        echo "PROCESS: SORT_BAM"
        echo "SAMPLE: ${sample}"
        date

        bash ${projectDir}/scripts/sort_bam.sh \
            ${bam} \
            ${sample}_sorted.bam \
            ${params.gatk}

    ) > ${logDir}/${sample}_sort.log 2>&1
    """
}


process MARK_DUPLICATES {
    tag "$sample"

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}_dedup.bam")

    script:
    """
    mkdir -p ${logDir}

    (
        echo "PROCESS: MARK_DUPLICATES"
        echo "SAMPLE: ${sample}"
        date

        bash ${projectDir}/scripts/mark_duplicates.sh \
            ${bam} \
            ${sample}_dedup.bam \
            ${sample}_metrics.txt \
            ${params.gatk}

    ) > ${logDir}/${sample}_markdup.log 2>&1
    """
}


process HAPLOTYPECALLER {
    tag "$sample"

    input:
    tuple val(sample), path(bam), val(fallback)

    output:
    tuple val(sample), path("${sample}.g.vcf.gz"), path("${sample}.g.vcf.gz.tbi")

    script:
    """
    mkdir -p ${logDir}

    (
        echo "PROCESS: HAPLOTYPECALLER"
        echo "SAMPLE: ${sample}"
        date

        ${params.gatk} --version

        bash ${projectDir}/scripts/haplotypecaller.sh \
            ${bam} \
            ${sample}.g.vcf.gz \
            ${params.ref} \
            ${params.gatk} \
            ${fallback}

    ) > ${logDir}/${sample}_haplotypecaller.log 2>&1
    """
}


process COMBINE_GVCFS {
    tag "combine"

    input:
    path gvcfs
    path indices

    output:
    path "combined.g.vcf.gz"
    path "combined.g.vcf.gz.tbi"

    script:
    def variants = gvcfs.collect { "--variant ${it}" }.join(" ")
    """
    mkdir -p ${logDir}

    (
        echo "PROCESS: COMBINE_GVCFS"
        date

        ${params.gatk} --version

        ${params.gatk} --java-options "-Xmx6g" CombineGVCFs \
            -R ${params.ref} \
            ${variants} \
            -O combined.g.vcf.gz

    ) > ${logDir}/combine_gvcfs.log 2>&1
    """
}


process GENOTYPE_GVCFS {
    tag "genotype"

    input:
    path combined
    path combined_tbi

    output:
    path "genotyped_variants.vcf.gz"
    path "genotyped_variants.vcf.gz.tbi"

    script:
    """
    mkdir -p ${logDir}

    (
        echo "PROCESS: GENOTYPE_GVCFS"
        date

        ${params.gatk} --version

        ${params.gatk} --java-options "-Xmx6g" GenotypeGVCFs \
            -R ${params.ref} \
            -V ${combined} \
            -O genotyped_variants.vcf.gz

    ) > ${logDir}/genotype.log 2>&1
    """
}


process FILTER_VARIANTS {
    tag "filter"

    input:
    path vcf
    path vcf_tbi

    output:
    path "filtered/filtered_snps.vcf", emit: snps
    path "filtered/filtered_indels.vcf", emit: indels

    script:
    """
    mkdir -p ${logDir}

    (
        echo "PROCESS: FILTER_VARIANTS"
        date

        bash ${projectDir}/scripts/filter_variants.sh \
            ${params.ref} \
            ${vcf} \
            filtered \
            ${params.gatk}

    ) > ${logDir}/filter.log 2>&1
    """
}


process ANNOTATE {
    tag "annotate"

    input:
    path snps_vcf
    path indels_vcf

    output:
    path "annotation/filtered_snps.hg38_multianno.txt", emit: snps_txt
    path "annotation/filtered_snps.hg38_multianno.vcf", emit: snps_vcf
    path "annotation/filtered_indels.hg38_multianno.txt", emit: indels_txt
    path "annotation/filtered_indels.hg38_multianno.vcf", emit: indels_vcf

    script:
    """
    mkdir -p ${logDir} annotation

    (
        echo "PROCESS: ANNOTATE"
        date

        perl ${params.annovar}/table_annovar.pl \
            ${snps_vcf} \
            ${params.annovar}/humandb \
            -buildver hg38 \
            -out annotation/filtered_snps \
            -remove \
            -protocol refGeneWithVer,avsnp150 \
            -operation g,f \
            -nastring . \
            --vcfinput \
            -thread 4 || true

        if [ ! -f annotation/filtered_snps.hg38_multianno.txt ]; then
            echo "ERROR: ANNOVAR failed on SNPs"
            exit 1
        fi

        perl ${params.annovar}/table_annovar.pl \
            ${indels_vcf} \
            ${params.annovar}/humandb \
            -buildver hg38 \
            -out annotation/filtered_indels \
            -remove \
            -protocol refGeneWithVer,avsnp150 \
            -operation g,f \
            -nastring . \
            --vcfinput \
            -thread 4 || true

        if [ ! -f annotation/filtered_indels.hg38_multianno.txt ]; then
            echo "ERROR: ANNOVAR failed on indels"
            exit 1
        fi

    ) > ${logDir}/annotate.log 2>&1
    """
}


process ADD_RSIDS {
    tag "rsids"
    publishDir "${params.results}/${workflow.runName}/annotation", mode: 'copy'

    input:
    path snps_vcf
    path indels_vcf

    output:
    path "filtered_snps.rsid.vcf"
    path "filtered_indels.rsid.vcf"

    script:
    """
    mkdir -p ${logDir}

    (
        echo "PROCESS: ADD_RSIDS"
        date

        bcftools annotate \
            --set-id +'%INFO/avsnp150' \
            ${snps_vcf} \
            -o filtered_snps.rsid.vcf

        bcftools annotate \
            --set-id +'%INFO/avsnp150' \
            ${indels_vcf} \
            -o filtered_indels.rsid.vcf

    ) > ${logDir}/add_rsids.log 2>&1
    """
}


process BEAUTIFY_TABLE {
    tag "beautify"
    publishDir "${params.results}/${workflow.runName}/summary", mode: 'copy'

    input:
    path snps_txt
    path indels_txt

    output:
    path "readable_snps.txt"
    path "readable_indels.txt"
    path "${params.gene}_report.txt"

    script:
    def sample1 = params.samples[0]
    def sample2 = params.samples[1]
    """
    mkdir -p ${logDir}

    (
        echo "PROCESS: BEAUTIFY_TABLE"
        date

        bash ${projectDir}/scripts/beautify_table.sh \
            ${snps_txt} \
            ${indels_txt} \
            ${params.gene} \
            ${sample1} \
            ${sample2}

    ) > ${logDir}/beautify.log 2>&1
    """
}


workflow {

    samples_ch = Channel.fromList(params.samples)

    bam_ch = samples_ch.map { sample ->
        tuple(sample, file("${params.project}/bam/${sample}.bam"))
    }

    fallback_ch = samples_ch.map { sample ->
        tuple(sample, "${params.fallback_dir}/${sample}.g.vcf.gz")
    }

    sorted = SORT_BAM(bam_ch)
    dedup  = MARK_DUPLICATES(sorted)

    dedup_with_fb = dedup.join(fallback_ch)

    gvcfs = HAPLOTYPECALLER(dedup_with_fb)

    all_gvcfs = gvcfs
        .map { sample, gvcf, tbi -> gvcf }
        .collect()

    all_tbis = gvcfs
        .map { sample, gvcf, tbi -> tbi }
        .collect()

    (combined, combined_tbi) = COMBINE_GVCFS(all_gvcfs, all_tbis)

    (genotype, genotype_tbi) = GENOTYPE_GVCFS(combined, combined_tbi)

    filtered = FILTER_VARIANTS(genotype, genotype_tbi)

    annotated = ANNOTATE(filtered.snps, filtered.indels)

    ADD_RSIDS(annotated.snps_vcf, annotated.indels_vcf)

    BEAUTIFY_TABLE(annotated.snps_txt, annotated.indels_txt)
}