nextflow.enable.dsl=2


/*
 * ============================================================
 * Parameters
 * ============================================================
 */

params.control_dir = "${baseDir}/data/control"
params.stim_dir    = "${baseDir}/data/stim"

params.outdir = "${baseDir}/results"

params.min_genes = 200
params.max_genes = 2500
params.max_mt = 5.0

params.expected_doublet_rate = 0.06

params.n_top_genes = 2000
params.n_pcs = 30


/*
 * ============================================================
 * Workflow
 * ============================================================
 */

workflow {


    /*
     * Input datasets
     *
     * Structure:
     *
     * data/
     * ├── control/
     * │   ├── matrix.mtx.gz
     * │   ├── barcodes.tsv.gz
     * │   └── features.tsv.gz
     * │
     * └── stim/
     *     ├── matrix.mtx.gz
     *     ├── barcodes.tsv.gz
     *     └── features.tsv.gz
     *
     */


    datasets = Channel.of(
        tuple(
            "control",
            file(params.control_dir)
        ),

        tuple(
            "stim",
            file(params.stim_dir)
        )
    )


    /*
     * QC filtering
     */

    qc_results = QC(
        datasets
    )


    /*
     * Doublet removal per condition
     */

    doublet_results = DOUBLETS(
        qc_results
    )


    /*
     * Merge conditions
     */

    merge_input = doublet_results
        .map { condition, file ->
            file
        }
        .collect()

    merged = MERGE(
        merge_input
    )


    /*
     * Normalization:
     *
     * normalize_total
     * log1p
     * highly_variable_genes
     * scale
     * PCA
     */

    normalized = NORMALIZE(
        merged
    )


    /*
     * Harmony integration
     */

    HARMONY(
        normalized
    )

}



/*
 * ============================================================
 * QC
 * ============================================================
 */


process QC {

    tag "${condition}"


    publishDir "${params.outdir}/qc/${condition}",
        mode: "copy",
        overwrite: true


    input:

    tuple val(condition), path(input_dir)


    output:

    tuple val(condition),
          path("qc_filtered.h5ad"),
          emit: qc


    script:


    """
    ${projectDir}/.venv/bin/python \
        ${projectDir}/scripts/qc.py \
        --input ${input_dir} \
        --output . \
        --condition ${condition} \
        --min-genes ${params.min_genes} \
        --max-genes ${params.max_genes} \
        --max-mt ${params.max_mt}
    """
}



/*
 * ============================================================
 * Doublet removal
 * ============================================================
 */


process DOUBLETS {

    tag "${condition}"


    publishDir "${params.outdir}/doublets/${condition}",
        mode: "copy",
        overwrite: true


    input:

    tuple val(condition), path(qc_file)


    output:
    tuple val(condition),
        path("${condition}_doublet_filtered.h5ad"),
        emit: doublets


    script:


    """
    ${projectDir}/.venv/bin/python \
        ${projectDir}/scripts/doublet_removal.py \
        --input ${qc_file} \
        --output . \
        --expected-doublet-rate ${params.expected_doublet_rate}

    mv doublet_filtered.h5ad ${condition}_doublet_filtered.h5ad
    """
}



/*
 * ============================================================
 * Merge
 * ============================================================
 */


process MERGE {

    tag "merge"


    publishDir "${params.outdir}/merged",
        mode: "copy",
        overwrite: true


    input:

    path(input_files)


    output:

    path("merged.h5ad"),
    emit: merged


    script:


    def inputs = input_files.join(" ")


   """
    ${projectDir}/.venv/bin/python \
        ${projectDir}/scripts/merge.py \
        --inputs ${inputs} \
        --output merged.h5ad
    """
}



/*
 * ============================================================
 * Normalize + PCA
 * ============================================================
 */


process NORMALIZE {

    tag "normalize"


    publishDir "${params.outdir}/normalized",
        mode: "copy",
        overwrite: true


    input:

    path(merged)


    output:

    path("normalized_pca.h5ad"),
    emit: normalized


    script:


    """
    ${projectDir}/.venv/bin/python \
        ${projectDir}/scripts/normalize.py \
        --input ${merged} \
        --output . \
        --n-top-genes ${params.n_top_genes} \
        --n-pcs ${params.n_pcs}
    """
}



/*
 * ============================================================
 * Harmony
 * ============================================================
 */


process HARMONY {

    tag "harmony"


    publishDir "${params.outdir}/integration/harmony",
        mode: "copy",
        overwrite: true


    input:

    path(normalized)


    output:

    path("harmony.h5ad"),
    emit: integrated


    script:


    """
    ${projectDir}/.venv/bin/python \
        ${projectDir}/scripts/harmony.py \
        --input ${normalized} \
        --output . \
        --n-pcs ${params.n_pcs}
    """
}