nextflow.enable.dsl = 2

/*
 * ---------------------------------------------------------------------
 *  Rooted maximum-likelihood phylogeny of the Caucasian grapevine set
 *
 *      413-sample rooted VCF  ->  preflight checks
 *                            ->  SNPhylo: LD pruning + SNP sequence
 *                            ->  PHYLIP dnaml, rooted on ZZ01
 *                            ->  phangorn bootstrap, 100 replicates
 *                            ->  validation of both trees
 *                            ->  tip annotation from the K=7 Q matrix
 *                            ->  fan and rectangular figures
 *
 *  Run:
 *      cd pipelines/05_tree
 *      nextflow run .
 *
 *      # reuse an ML tree that was produced outside the pipeline
 *      nextflow run . --ml_tree /path/to/cauc_rooted.ml.tree \
 *                     --phylip  /path/to/cauc_rooted.phylip.txt
 *
 *  The ML stage is one long single-threaded dnaml run; the bootstrap is
 *  100 independent replicates. They are separate processes so that the
 *  second can ask for many cores that the first would only leave idle.
 * ---------------------------------------------------------------------
 */


params.tree_outdir = "${params.outdir}/tree"


process PREFLIGHT {

    tag "413 samples, ZZ01, ${params.expect_variants} variants"


    publishDir "${params.tree_outdir}",
        mode: "copy"


    input:

    path vcf_gz


    output:

    path "preflight.txt"


    script:

    """

    preflight_vcf.py \
        --vcf ${vcf_gz} \
        --expect-samples ${params.expect_samples} \
        --expect-variants ${params.expect_variants} \
        --outgroup ${params.outgroup} \
        --out preflight.txt

    """

}


process SNPHYLO_ML {

    tag "dnaml, outgroup ${params.outgroup}"


    publishDir "${params.tree_outdir}/snphylo",
        mode: "copy"


    input:

    path vcf_gz

    path preflight


    output:

    path "${params.prefix}.ml.tree",    emit: ml_tree

    path "${params.prefix}.phylip.txt", emit: phylip

    path "${params.prefix}.ml.txt"

    path "${params.prefix}.fasta"

    path "${params.prefix}.id.txt"

    path "run_summary.txt"


    script:

    //
    // SNPhylo insists on an uncompressed VCF. It is decompressed into
    // the task's own working directory -- never next to the shared
    // input, which the task brief forbids.
    //
    // snphylo.vcf.sh is copied in rather than run from the shared tree.
    // bash reads a script incrementally as it executes it, and this one
    // is open for the six hours dnaml takes; on the project NFS mount
    // the handle went stale partway through and bash died with
    //
    //     snphylo.vcf.sh: error reading input file: Stale file handle
    //
    // seconds after dnaml had finished, discarding the tree. A local
    // copy is not on NFS and cannot go stale.
    //

    """

    cp ${params.snphylo} ./snphylo.sh

    gunzip -c ${vcf_gz} > input.vcf

    status=0

    bash ./snphylo.sh \
        -v input.vcf \
        -r \
        -m ${params.maf} \
        -M ${params.missing} \
        -l ${params.ld} \
        -a ${params.last_autosome} \
        -o ${params.outgroup} \
        -t ${task.cpus} \
        -P ${params.prefix} || status=\$?

    #
    # dnaml writes "outfile" and "outtree"; snphylo.vcf.sh renames them
    # as its last act, after an R plotting step that is allowed to fail.
    # Recover them here so a cosmetic failure at the end does not throw
    # away six hours of inference.
    #

    if [ ! -s ${params.prefix}.ml.tree ] && [ -s outtree ]; then
        echo "recovering the tree from PHYLIP's outtree" >&2
        mv outtree ${params.prefix}.ml.tree
    fi

    if [ ! -s ${params.prefix}.ml.txt ] && [ -s outfile ]; then
        mv outfile ${params.prefix}.ml.txt
    fi

    if [ ! -s ${params.prefix}.ml.tree ]; then
        echo "snphylo exited \$status and left no tree" >&2
        exit 1
    fi

    rm -f input.vcf infile snphylo.sh

    {
        echo "input           ${params.vcf_rooted}"
        echo "snphylo         ${params.snphylo}"
        echo "maf             >= ${params.maf}"
        echo "missing rate    <= ${params.missing}"
        echo "ld threshold    ${params.ld}"
        echo "chromosomes     1-${params.last_autosome}"
        echo "outgroup        ${params.outgroup}"
        echo "low-depth screen skipped (-r)"
        echo "sites retained  \$(head -1 ${params.prefix}.phylip.txt | awk '{print \$2}')"
        echo "slurm job       \${SLURM_JOB_ID:-none}"
        echo "snphylo exit    \$status"
    } > run_summary.txt

    """

}


process BOOTSTRAP {

    tag "${params.bootstrap} replicates on ${task.cpus} cores"


    publishDir "${params.tree_outdir}/snphylo",
        mode: "copy"


    input:

    path ml_tree

    path phylip


    output:

    path "${params.prefix}.bs.tree", emit: bs_tree


    script:

    //
    // SNPhylo's own determine_bs_tree.R calls bootstrap.pml() with
    // multicore=TRUE and no mc.cores, so phangorn falls back to
    // detectCores() -- 64 on these nodes, whatever SLURM actually
    // allocated. That oversubscription is what killed the group's
    // earlier attempt. Here the count is passed in explicitly.
    //

    """

    export R_LIBS=${params.snphylo_r_libs}

    bootstrap_tree.R \
        --tree ${ml_tree} \
        --phylip ${phylip} \
        --replicates ${params.bootstrap} \
        --cores ${task.cpus} \
        --out ${params.prefix}.bs.tree

    """

}


process VALIDATE {

    tag "413 tips, ZZ01 present"


    publishDir "${params.tree_outdir}",
        mode: "copy"


    input:

    path ml_tree

    path bs_tree


    output:

    path "tree_validation.txt"


    script:

    """

    validate_tree.py \
        --ml-tree ${ml_tree} \
        --bs-tree ${bs_tree} \
        --expect-tips ${params.expect_samples} \
        --outgroup ${params.outgroup} \
        --out tree_validation.txt

    """

}


process ANNOTATE {

    tag "K=${params.k}, Q >= ${params.min_q}"


    publishDir "${params.tree_outdir}",
        mode: "copy"


    input:

    path bs_tree

    path q_file

    path fam

    path metadata


    output:

    path "tip_annotation.tsv", emit: annotation

    path "group_labels.tsv",   emit: labels


    script:

    """

    annotate_tips.py \
        --tree ${bs_tree} \
        --q ${q_file} \
        --fam ${fam} \
        --metadata ${metadata} \
        --min-q ${params.min_q} \
        --outgroup ${params.outgroup} \
        --outdir .

    """

}


process PLOT {

    tag "${layout}"


    publishDir "${params.tree_outdir}/plots",
        mode: "copy"


    input:

    val layout

    path bs_tree

    path annotation

    path labels


    output:

    path "tree_*.pdf"

    path "tree_*.png"


    script:

    """

    plot_tree.py \
        --tree ${bs_tree} \
        --annotation ${annotation} \
        --labels ${labels} \
        --layout ${layout} \
        --outgroup ${params.outgroup} \
        --min-support ${params.min_support} \
        --outdir .

    """

}


workflow {

    vcf_gz = Channel.value(
        file(params.vcf_rooted, checkIfExists: true)
    )

    q_file = Channel.value(
        file(params.q_file, checkIfExists: true)
    )

    fam = Channel.value(
        file(params.fam, checkIfExists: true)
    )

    metadata = Channel.value(
        file(params.metadata, checkIfExists: true)
    )


    PREFLIGHT(vcf_gz)


    /*
     * dnaml on 413 taxa takes hours. When it has already been run --
     * by hand, or by an earlier attempt of this pipeline -- point
     * --ml_tree and --phylip at its output and start from the
     * bootstrap.
     */

    if (params.ml_tree || params.phylip) {

        if (!params.ml_tree || !params.phylip) {

            error """
            --ml_tree and --phylip go together.

            The bootstrap needs both the starting topology and the
            alignment it was inferred from; a tree resampled against a
            different alignment is wrong in a way nothing downstream
            would catch.

              --ml_tree ${params.ml_tree ?: '(not given)'}
              --phylip  ${params.phylip ?: '(not given)'}
            """.stripIndent()
        }

        ml_tree = Channel.value(file(params.ml_tree, checkIfExists: true))

        phylip = Channel.value(file(params.phylip, checkIfExists: true))

    }
    else {

        SNPHYLO_ML(vcf_gz, PREFLIGHT.out)

        ml_tree = SNPHYLO_ML.out.ml_tree

        phylip = SNPHYLO_ML.out.phylip

    }


    BOOTSTRAP(ml_tree, phylip)

    VALIDATE(ml_tree, BOOTSTRAP.out.bs_tree)

    ANNOTATE(BOOTSTRAP.out.bs_tree, q_file, fam, metadata)

    layouts = Channel.fromList(
        params.layouts.toString().split(",").collect { it.trim() }
    )

    PLOT(
        layouts,
        BOOTSTRAP.out.bs_tree,
        ANNOTATE.out.annotation,
        ANNOTATE.out.labels
    )

}


workflow.onComplete {

    log.info """
    ----------------------------------------------------------
    Rooted ML phylogeny
    input    : ${params.vcf_rooted}
    filters  : MAF >= ${params.maf}, missing <= ${params.missing}, LD ${params.ld}
    outgroup : ${params.outgroup}
    bootstrap: ${params.bootstrap} replicates
    results  : ${params.tree_outdir}
    status   : ${workflow.success ? 'OK' : 'FAILED'}
    ----------------------------------------------------------
    """.stripIndent()
}
