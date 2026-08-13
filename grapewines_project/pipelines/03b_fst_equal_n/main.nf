nextflow.enable.dsl = 2

/*
 * ---------------------------------------------------------------------
 *  Pairwise FST with every ancestry group cut to the same sample size
 *
 *      sample_q_values_K<K>.tsv  ->  R independent draws of n per group
 *                                ->  VCF cut to the assigned samples
 *                                ->  21 x R pairwise FST
 *                                ->  within-group pi, full and equal-n
 *                                ->  comparison against 03_fst
 *
 *  Why this exists
 *  ---------------
 *  In 03_fst the 21 FST values correlate with group size at Spearman
 *  -0.79. FST is between-group variance over total variance, and the
 *  groups are not comparable objects: at Q >= 0.75 a component with 12
 *  members contributes its tight core while one with 49 contributes a
 *  broad slice, so the denominator differs systematically with size.
 *
 *  This run holds the sample size fixed and asks whether the effect
 *  survives. It is a diagnostic, not a correction -- see the comment at
 *  the top of bin/draw_subsamples.py.
 *
 *  03_fst is left untouched: it is the reference run, and its results
 *  are what the group's figures are built on.
 *
 *  Run:
 *      cd pipelines/03b_fst_equal_n
 *      nextflow run .
 * ---------------------------------------------------------------------
 */


params.assignments =
    "${params.outdir}/fst/K${params.k}/sample_lists/sample_q_values_K${params.k}.tsv"

params.full_summary =
    "${params.outdir}/fst/K${params.k}/pairwise_fst/pairwise_fst_summary_K${params.k}.tsv"

params.group_sizes =
    "${params.outdir}/fst/K${params.k}/sample_lists/group_sizes_K${params.k}.tsv"

params.equal_n_outdir = "${params.outdir}/fst/K${params.k}_equal_n"


process DRAW_SUBSAMPLES {

    tag "n = ${params.subsample_n ?: 'smallest'}, ${params.replicates} draws"

    publishDir "${params.equal_n_outdir}/sample_lists",
        mode: "copy"


    input:

    path assignments


    output:

    path "rep*_K*_samples.txt", emit: groups

    path "assigned_samples.txt", emit: assigned


    script:

    """

    draw_subsamples.py \
        --assignments ${assignments} \
        --n ${params.subsample_n} \
        --replicates ${params.replicates} \
        --seed ${params.seed} \
        --outdir .

    """

}


/*
 * The full VCF is 2.1 GB and holds all 412 accessions. Every one of the
 * 21 x R vcftools passes would stream all of it. Cutting it once to the
 * 160 assigned samples makes each of those passes read a much smaller
 * file, which is the difference between a few hours and most of a day.
 */

process SUBSET_VCF {

    tag "160 assigned samples"


    input:

    path vcf

    path assigned


    output:

    path "assigned.vcf.gz"


    script:

    """

    bcftools view \
        --samples-file ${assigned} \
        --force-samples \
        --output-type z \
        --output assigned.vcf.gz \
        ${vcf}

    """

}


process PAIRWISE_FST {

    tag "${replicate}: ${pop_a.name.tokenize('_')[1]} vs ${pop_b.name.tokenize('_')[1]}"

    publishDir "${params.equal_n_outdir}/pairwise_fst",
        mode: "copy",
        pattern: "*.log"


    input:

    tuple val(replicate), path(pop_a), path(pop_b)

    path vcf


    output:

    path "*.log", emit: log


    script:

    def group_a = pop_a.name.replaceAll(/^rep\d+_/, "")
                            .replace("_samples.txt", "")

    def group_b = pop_b.name.replaceAll(/^rep\d+_/, "")
                            .replace("_samples.txt", "")

    def comparison = "${replicate}_${group_a}_vs_${group_b}"

    //
    // Only the log is kept. The per-SNP table is 800 000 rows and there
    // would be 210 of them; the mean is all this analysis needs, and
    // summarize_fst.py already established that vcftools' own mean and
    // the mean computed from the table agree.
    //
    // vcftools 0.1.13 prints its summary to stderr and writes no .log,
    // so the redirect below is what creates the file at all -- the same
    // fix as in 03_fst.
    //

    """

    vcftools \
        --gzvcf ${vcf} \
        --weir-fst-pop ${pop_a} \
        --weir-fst-pop ${pop_b} \
        --out ${comparison} \
        2> ${comparison}.stderr.txt

    cat ${comparison}.stderr.txt >&2

    if [ ! -s ${comparison}.log ]; then
        cp ${comparison}.stderr.txt ${comparison}.log
    fi

    #
    # 210 of these tables at ~30 MB each would be 6 GB of work
    # directory on a 30 GB NFS home, and none of it is read again --
    # only the mean in the log is. Drop it as soon as the log exists.
    #
    rm -f ${comparison}.weir.fst

    """

}


/*
 * Nucleotide diversity inside each group. This is the denominator of
 * FST measured directly, so it says whether the small groups are narrow
 * because of how few samples they have or because of what they are.
 */

process GROUP_PI {

    tag "${scope} ${group}"

    publishDir "${params.equal_n_outdir}/within_group_pi",
        mode: "copy"


    input:

    tuple val(scope), val(group), path(samples)

    path vcf


    output:

    path "${scope}_${group}.pi"


    script:

    """

    vcftools \
        --gzvcf ${vcf} \
        --keep ${samples} \
        --site-pi \
        --out ${scope}_${group} \
        2>&1

    awk 'NR > 1 && \$3 != "-nan" { total += \$3; n++ }
         END { if (n > 0) printf "%.8f\\n", total / n; else print "" }' \
        ${scope}_${group}.sites.pi > ${scope}_${group}.pi

    if [ ! -s ${scope}_${group}.pi ]; then
        echo "no usable per-site pi for ${scope} ${group}" >&2
        exit 1
    fi

    rm -f ${scope}_${group}.sites.pi

    """

}


process COLLECT {

    publishDir "${params.equal_n_outdir}",
        mode: "copy"


    input:

    path logs

    path pis

    path full_summary

    path group_sizes


    output:

    path "*.tsv"


    script:

    """

    mkdir -p fst pi

    for f in ${logs}; do cp "\$f" fst/; done
    for f in ${pis};  do cp "\$f" pi/;  done

    collect_equal_n.py \
        --fstdir fst \
        --pidir pi \
        --full-summary ${full_summary} \
        --group-sizes ${group_sizes} \
        --outdir .

    """

}


workflow {

    log.info """
    ----------------------------------------------------------
    FST with equal group sizes, K = ${params.k}
    draws   : ${params.replicates} x ${params.subsample_n ?: 'smallest group'} per group
    results : ${params.equal_n_outdir}
    ----------------------------------------------------------
    """.stripIndent()


    assignments = Channel.fromPath(params.assignments, checkIfExists: true)

    vcf = Channel.fromPath(params.vcf, checkIfExists: true)


    DRAW_SUBSAMPLES(assignments)

    SUBSET_VCF(vcf, DRAW_SUBSAMPLES.out.assigned)

    small_vcf = SUBSET_VCF.out.first()


    /*
     * The draw process emits every list at once. Regroup by replicate,
     * then form the 21 unordered pairs inside each replicate.
     */

    by_replicate = DRAW_SUBSAMPLES.out.groups
        .flatten()
        .map { file ->
            tuple((file.name =~ /^(rep\d+)_/)[0][1], file)
        }
        .groupTuple()


    pairs = by_replicate
        .flatMap { replicate, files ->

            def sorted = files.sort { it.name }

            def combinations = []

            for (int i = 0; i < sorted.size() - 1; i++) {
                for (int j = i + 1; j < sorted.size(); j++) {
                    combinations << tuple(replicate, sorted[i], sorted[j])
                }
            }

            combinations
        }


    PAIRWISE_FST(pairs, small_vcf)


    /*
     * pi is measured twice: on the full groups, which is the honest
     * description of each group, and on every equal-n draw, which says
     * whether the differences are only about sample size.
     */

    full_groups = Channel
        .fromPath(
            "${params.outdir}/fst/K${params.k}/sample_lists/K*_samples.txt",
            checkIfExists: true
        )
        .map { file ->
            tuple("full", file.name.replace("_samples.txt", ""), file)
        }

    drawn_groups = DRAW_SUBSAMPLES.out.groups
        .flatten()
        .map { file ->
            def matcher = (file.name =~ /^(rep\d+)_(K\d+)_samples\.txt$/)
            tuple(matcher[0][1], matcher[0][2], file)
        }

    GROUP_PI(full_groups.mix(drawn_groups), small_vcf)


    COLLECT(
        PAIRWISE_FST.out.log.collect(),
        GROUP_PI.out.collect(),
        Channel.fromPath(params.full_summary, checkIfExists: true),
        Channel.fromPath(params.group_sizes, checkIfExists: true)
    )
}


workflow.onComplete {

    log.info """
    ----------------------------------------------------------
    FST with equal group sizes, K = ${params.k}
    results : ${params.equal_n_outdir}
    status  : ${workflow.success ? 'OK' : 'FAILED'}
    ----------------------------------------------------------
    """.stripIndent()
}
