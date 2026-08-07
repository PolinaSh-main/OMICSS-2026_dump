nextflow.enable.dsl=2


params.results_dir =
    "/mnt/nas0/user/polina.shevyakova/grapewines_project/admixture/results"

params.scripts_dir =
    "/mnt/nas0/user/polina.shevyakova/grapewines_project/admixture/scripts"

params.metadata =
    "/mnt/nas1/proj/omicss26/gp3/data/metadata/cauc_grape_metadata.csv"

params.fam =
    "/mnt/nas0/user/polina.shevyakova/grapewines_project/admixture/data_filtered/cauc_filtered.final.fam"

params.output_dir =
    "/mnt/nas0/user/polina.shevyakova/grapewines_project/admixture/output"

params.admixture_log_dir =
    "/mnt/nas0/user/polina.shevyakova/grapewines_project/admixture/results/log"



process EXTRACT_CV {

    tag "Extract ADMIXTURE CV"


    publishDir "${params.output_dir}/cv",
        mode: "copy"


    input:

    path logdir


    output:

    path "cv.tsv"


    script:

    """

    python3 ${params.scripts_dir}/extract_cv.py \
        --logdir ${logdir} \
        --output cv.tsv

    """

}



process PLOT_CV {

    tag "Plot CV curve"


    publishDir "${params.output_dir}/cv",
        mode: "copy"


    input:

    path cv


    output:

    path "cv_plot.png"

    path "cv_plot.pdf"



    script:

    """

    python3 ${params.scripts_dir}/plot_cv.py \
        --input ${cv} \
        --png cv_plot.png \
        --pdf cv_plot.pdf

    """

}



process BUILD_ORDERS {


    tag "Build ADMIXTURE orders"


    publishDir "${params.output_dir}/orders",
        mode: "copy"



    input:

    path q_files

    path fam

    path metadata



    output:

    path "order_Q*.tsv"



    script:

    """

    python3 ${params.scripts_dir}/build_orders.py \
        --qdir . \
        --fam ${fam} \
        --metadata ${metadata} \
        --outdir .

    """

}



process PLOT_ADMIXTURE_GRID {


    tag "ADMIXTURE ordering comparison"


    publishDir "${params.output_dir}/admixture",
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

    python3 ${params.scripts_dir}/plot_admixture_grid.py \
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
     * ADMIXTURE Q files  →  value channel via .collect()
     */

    q_collected = Channel
        .fromPath(
            "${params.results_dir}/cauc_filtered.final.*.Q",
            checkIfExists: true
        )
        .collect()



    /*
     * Single-file inputs  →  value channels (reusable)
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