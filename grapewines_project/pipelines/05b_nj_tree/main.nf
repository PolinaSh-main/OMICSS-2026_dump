nextflow.enable.dsl = 2

/*
 * ---------------------------------------------------------------------
 *  Neighbour-joining tree from the filtered SNP set
 *
 *      PLINK bed/bim/fam  ->  LD pruning
 *                         ->  1-IBS distance + genotype export
 *                         ->  NJ tree, bootstrapped, midpoint rooted
 *                         ->  circular tree coloured by metadata
 *
 *  Run:
 *      cd pipelines/05_tree
 *      nextflow run . --bfile /path/to/cauc_filtered.final
 *
 *      # colour by ADMIXTURE group instead of country
 *      nextflow run . --color_by assignment \
 *          --groups ../../results/fst/K7/sample_lists/sample_q_values_K7.tsv
 *
 *  LD pruning first is not optional here. Linked SNPs carry the same
 *  information several times over, which distorts branch lengths and
 *  makes bootstrap support look far better than it is.
 * ---------------------------------------------------------------------
 */


params.tree_outdir = "${params.outdir}/tree"


process LD_PRUNE {

    tag "r2 < ${params.ld_r2}"


    publishDir "${params.tree_outdir}",
        mode: "copy",
        pattern: "*.log"


    input:

    tuple path(bed), path(bim), path(fam)


    output:

    path "pruned.prune.in", emit: keep

    path "*.log"


    script:

    def prefix = bed.baseName

    """

    plink \
        --bfile ${prefix} \
        --allow-extra-chr 0 \
        --chr-set ${params.chr_set} \
        --indep-pairwise ${params.ld_window} ${params.ld_step} ${params.ld_r2} \
        --out pruned

    """

}


process DISTANCE_MATRIX {

    tag "1 - IBS"


    publishDir "${params.tree_outdir}",
        mode: "copy",
        pattern: "*.log"


    input:

    tuple path(bed), path(bim), path(fam)

    path keep


    output:

    path "tree_input.traw",      emit: traw

    path "tree_input.mdist",     emit: mdist

    path "tree_input.mdist.id",  emit: ids

    path "*.log"


    script:

    def prefix = bed.baseName

    //
    // Two calls rather than one: PLINK 1.9 does not accept --distance
    // and --recode in the same invocation.
    //

    """

    plink \
        --bfile ${prefix} \
        --allow-extra-chr 0 \
        --chr-set ${params.chr_set} \
        --extract ${keep} \
        --distance square 1-ibs \
        --out tree_input

    plink \
        --bfile ${prefix} \
        --allow-extra-chr 0 \
        --chr-set ${params.chr_set} \
        --extract ${keep} \
        --recode A-transpose \
        --out tree_input_geno

    mv tree_input_geno.traw tree_input.traw

    """

}


process BUILD_TREE {

    tag "NJ, ${params.bootstrap} bootstrap replicates"


    publishDir "${params.tree_outdir}",
        mode: "copy"


    input:

    path traw

    path mdist


    output:

    path "tree.nwk", emit: tree

    path "distance_matrix.tsv"

    path "tree_tips.tsv"


    script:

    def outgroup = params.outgroup ? "--outgroup ${params.outgroup}" : ""

    """

    build_tree.py \
        --traw ${traw} \
        --mdist ${mdist} \
        --bootstrap ${params.bootstrap} \
        --root ${params.root} \
        ${outgroup} \
        --seed ${params.seed} \
        --outdir .

    """

}


process PLOT_TREE {

    tag "${colour_by} / ${scale}"


    publishDir "${params.tree_outdir}/plots",
        mode: "copy"


    input:

    tuple val(colour_by), val(scale)

    path tree

    path metadata

    path groups


    output:

    path "tree_*.pdf"

    path "tree_*.png"


    script:

    def labels = params.show_labels ? "--show-labels" : ""

    def group_file = groups.name != "NO_GROUPS" ? "--groups ${groups}" : ""

    """

    plot_tree.py \
        --tree ${tree} \
        --metadata ${metadata} \
        ${group_file} \
        --color-by '${colour_by}' \
        --layout ${params.layout} \
        --scale ${scale} \
        ${labels} \
        --outdir .

    """

}


workflow {

    /*
     * PLINK wants all three files staged next to each other.
     */

    plink_set = Channel.value(
        [
            file("${params.bfile}.bed", checkIfExists: true),
            file("${params.bfile}.bim", checkIfExists: true),
            file("${params.bfile}.fam", checkIfExists: true),
        ]
    )

    metadata = Channel.value(
        file(params.metadata, checkIfExists: true)
    )

    /*
     * Optional: colouring by ADMIXTURE group needs 03_fst output. A
     * placeholder keeps the process signature the same when it is not
     * wanted.
     */

    groups = Channel.value(
        params.groups
            ? file(params.groups, checkIfExists: true)
            : file("${projectDir}/assets/NO_GROUPS")
    )


    LD_PRUNE(plink_set)

    DISTANCE_MATRIX(plink_set, LD_PRUNE.out.keep)

    BUILD_TREE(
        DISTANCE_MATRIX.out.traw,
        DISTANCE_MATRIX.out.mdist
    )

    /*
     * One figure per colouring per scale, in parallel.
     *
     * Both scales are drawn by default: the linear one is the honest
     * picture of the distances, the cladogram is the one the topology
     * is actually readable on.
     */

    colourings = Channel.fromList(
        params.color_by.toString().split(",").collect { it.trim() }
    )

    scales = Channel.fromList(
        params.scale.toString().split(",").collect { it.trim() }
    )

    PLOT_TREE(
        colourings.combine(scales),
        BUILD_TREE.out.tree,
        metadata,
        groups
    )

}


workflow.onComplete {

    log.info """
    ----------------------------------------------------------
    Neighbour-joining tree
    distance : 1 - IBS on LD-pruned SNPs (r2 < ${params.ld_r2})
    bootstrap: ${params.bootstrap} replicates
    rooting  : ${params.root}
    results  : ${params.tree_outdir}
    status   : ${workflow.success ? 'OK' : 'FAILED'}
    ----------------------------------------------------------
    """.stripIndent()
}
