nextflow.enable.dsl = 2

/*
 * ---------------------------------------------------------------------
 *  ADMIXTURE post-processing
 *
 *  ADMIXTURE itself is still run by bin/run_admixture.sh (a plain
 *  SLURM loop over K). This pipeline takes its output and produces
 *  everything that goes into the report:
 *
 *      *.Q + logs  ->  CV table  ->  CV curve
 *                  ->  sample orderings  ->  barplot grid
 *
 *  Helper scripts live in bin/ and are on PATH automatically.
 * ---------------------------------------------------------------------
 */


process EXTRACT_CV {

    tag "Extract ADMIXTURE CV"


    publishDir "${params.outdir}/cv",
        mode: "copy"


    input:

    path logdir


    output:

    path "cv.tsv"


    script:

    """

    extract_cv.py \
        --logdir ${logdir} \
        --output cv.tsv

    """

}



process PLOT_CV {

    tag "Plot CV curve"


    publishDir "${params.outdir}/cv",
        mode: "copy"


    input:

    path cv


    output:

    path "cv_plot.png"

    path "cv_plot.pdf"



    script:

    """

    plot_cv.py \
        --input ${cv} \
        --png cv_plot.png \
        --pdf cv_plot.pdf

    """

}



process BUILD_ORDERS {


    tag "Build ADMIXTURE orders"


    publishDir "${params.outdir}/orders",
        mode: "copy"



    input:

    path q_files

    path fam

    path metadata



    output:

    path "order_Q*.tsv"



    script:

    """

    build_orders.py \
        --qdir . \
        --fam ${fam} \
        --metadata ${metadata} \
        --outdir .

    """

}



process PLOT_ADMIXTURE_GRID {


    tag "ADMIXTURE ordering comparison"


    publishDir "${params.outdir}/orderings",
        mode: "copy"



    input:

    path orders

    path q_files

    path fam

    path metadata



    output:

    path "*.png"



    script:

    """

    plot_admixture_grid.py \
        --qdir . \
        --orders . \
        --fam ${fam} \
        --metadata ${metadata} \
        --outdir .

    """

}



workflow {


    /*
     * CV extraction
     */

    logdir = Channel
        .fromPath(
            params.admixture_log_dir,
            checkIfExists: true
        )


    cv = EXTRACT_CV(
        logdir
    )


    PLOT_CV(
        cv
    )



    /*
     * ADMIXTURE Q files  ->  value channel via .collect()
     */

    q_collected = Channel
        .fromPath(
            "${params.admixture_dir}/cauc_filtered.final.*.Q",
            checkIfExists: true
        )
        .collect()



    /*
     * Single-file inputs  ->  value channels (reusable)
     */

    fam = Channel.value(
        file(params.fam, checkIfExists: true)
    )


    metadata = Channel.value(
        file(params.metadata, checkIfExists: true)
    )



    /*
     * Build sample orders:
     *
     * order_Q2.tsv
     * ...
     * order_Q10.tsv
     */

    orders = BUILD_ORDERS(
        q_collected,
        fam,
        metadata
    )



    /*
     * Plot all orderings
     */

    PLOT_ADMIXTURE_GRID(
        orders,
        q_collected,
        fam,
        metadata
    )

}
