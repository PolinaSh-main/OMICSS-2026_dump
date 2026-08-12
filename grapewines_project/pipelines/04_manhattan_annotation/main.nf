nextflow.enable.dsl = 2

/*
 * ---------------------------------------------------------------------
 *  From per-SNP FST to candidate genes
 *
 *      *.weir.fst  ->  50 kb windows, mean FST per window
 *                  ->  Manhattan plot, top 1% windows marked
 *                  ->  adjacent outlier windows merged into regions
 *                  ->  genes in and around those regions
 *                  ->  a short, ranked table for the report
 *
 *  Input is whatever pipelines/03_fst already produced, so nothing is
 *  recomputed from the VCF.
 *
 *  Run:
 *      cd pipelines/04_manhattan_annotation
 *      nextflow run . --k 7
 *
 *      # a single comparison
 *      nextflow run . --k 7 --comparisons K1_vs_K3
 * ---------------------------------------------------------------------
 */


params.fst_dir = "${params.outdir}/fst/K${params.k}/pairwise_fst"

params.manhattan_outdir = "${params.outdir}/manhattan/K${params.k}"


process WINDOW_FST {

    tag "${comparison}"


    publishDir "${params.manhattan_outdir}/${comparison}",
        mode: "copy"


    input:

    tuple val(comparison), path(fst)


    output:

    tuple val(comparison), path("*_outlier_windows.tsv"), emit: outliers

    path "*_windows.tsv"

    path "*_manhattan.pdf"

    path "*_manhattan.png"


    script:

    """

    window_fst.py \
        --fst ${fst} \
        --comparison ${comparison} \
        --window ${params.window} \
        --min-snps ${params.min_snps} \
        --top-fraction ${params.top_fraction} \
        --outdir .

    """

}


process MERGE_REGIONS {

    tag "${comparison}"


    publishDir "${params.manhattan_outdir}/${comparison}",
        mode: "copy"


    input:

    tuple val(comparison), path(outliers)


    output:

    tuple val(comparison), path("*_candidate_regions.tsv")


    script:

    """

    merge_regions.py \
        --windows ${outliers} \
        --comparison ${comparison} \
        --max-gap ${params.max_gap} \
        --outdir .

    """

}


process ANNOTATE_REGIONS {

    tag "${comparison}"


    publishDir "${params.manhattan_outdir}/${comparison}",
        mode: "copy"


    input:

    tuple val(comparison), path(regions)

    path gff

    path b2g


    output:

    tuple val(comparison), path("*_candidate_genes.tsv")


    script:

    """

    annotate_regions.py \
        --regions ${regions} \
        --gff ${gff} \
        --b2g ${b2g} \
        --flank ${params.flank} \
        --comparison ${comparison} \
        --outdir .

    """

}


process SELECT_GENES {

    tag "${comparison}"


    publishDir "${params.manhattan_outdir}/${comparison}",
        mode: "copy"


    input:

    tuple val(comparison), path(genes)


    output:

    path "*_selected_regions.tsv"

    path "*_low_information_genes.tsv"


    script:

    """

    select_candidate_genes.py \
        --genes ${genes} \
        --top-regions ${params.top_regions} \
        --genes-per-region ${params.genes_per_region} \
        --comparison ${comparison} \
        --outdir .

    """

}


workflow {

    /*
     * "all", or a comma-separated list such as "K1_vs_K3,K2_vs_K5".
     */

    def wanted = params.comparisons.toString().toLowerCase() == "all"
        ? null
        : params.comparisons.toString().split(",").collect { it.trim() }.toSet()


    fst_files = Channel
        .fromPath(
            "${params.fst_dir}/*.weir.fst",
            checkIfExists: true
        )
        .map { path ->

            def comparison = path.name.replace(".weir.fst", "")

            return [ comparison, path ]
        }
        .filter { item ->

            wanted == null || wanted.contains(item[0])
        }


    gff = Channel.value(
        file(params.gff, checkIfExists: true)
    )

    b2g = Channel.value(
        file(params.b2g, checkIfExists: true)
    )


    WINDOW_FST(fst_files)

    MERGE_REGIONS(WINDOW_FST.out.outliers)

    ANNOTATE_REGIONS(
        MERGE_REGIONS.out,
        gff,
        b2g
    )

    SELECT_GENES(ANNOTATE_REGIONS.out)

}


workflow.onComplete {

    log.info """
    ----------------------------------------------------------
    Manhattan + gene annotation, K = ${params.k}
    windows : ${params.window} bp, min ${params.min_snps} SNPs
    outliers: top ${params.top_fraction * 100}%
    results : ${params.manhattan_outdir}
    status  : ${workflow.success ? 'OK' : 'FAILED'}
    ----------------------------------------------------------
    """.stripIndent()
}
