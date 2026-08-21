nextflow.enable.dsl = 2

/*
 * ---------------------------------------------------------------------
 *  Population differentiation between ADMIXTURE groups (FST)
 *
 *      .Q + .fam  ->  sample lists per ancestry group
 *                 ->  metadata profile + proposed group names
 *                 ->  pairwise FST for every pair of groups
 *                 ->  summary table
 *                 ->  heatmap
 *
 *  Run:
 *      cd pipelines/03_fst
 *      nextflow run . --k 7
 *
 *      # one comparison only
 *      nextflow run . --k 7 --pairs K1:K5
 *
 *  At K=7 this is 21 comparisons. They are independent, so each one is
 *  its own SLURM job and they run in parallel rather than in the serial
 *  loop the reference scripts use.
 * ---------------------------------------------------------------------
 */


/*
 * Derived from --k, so it has to be built here rather than in
 * nextflow.config: values given on the command line are applied after
 * the config is parsed, and a config-time GString would have baked in
 * the default K instead.
 */

params.q_file = "${params.admixture_dir}/cauc_filtered.final.${params.k}.Q"

params.fst_outdir = "${params.outdir}/fst/K${params.k}"


process ASSIGN_GROUPS {

    tag "K=${params.k}, min Q=${params.min_q}"


    publishDir "${params.fst_outdir}/sample_lists",
        mode: "copy"


    input:

    path fam

    path q_file


    output:

    path "K*_samples.txt",         emit: lists

    path "sample_q_values_K*.tsv", emit: q_values

    path "group_sizes_K*.tsv",     emit: sizes

    path "admixed_samples.txt"


    script:

    """

    assign_admixture_groups.py \
        --fam ${fam} \
        --q ${q_file} \
        --min-q ${params.min_q} \
        --outdir .

    """

}


process DESCRIBE_GROUPS {

    tag "Name K1..K${params.k} from metadata"


    publishDir "${params.fst_outdir}/sample_lists",
        mode: "copy"


    input:

    path q_values

    path metadata


    output:

    path "group_labels_K*.tsv", emit: labels

    path "group_interpretation_K*.tsv"

    path "group_composition_K*.tsv"


    script:

    """

    describe_groups.py \
        --assignments ${q_values} \
        --metadata ${metadata} \
        --outdir .

    """

}


process PAIRWISE_FST {

    tag "${pop_a.simpleName.replace('_samples', '')} vs ${pop_b.simpleName.replace('_samples', '')}"


    publishDir "${params.fst_outdir}/pairwise_fst",
        mode: "copy"


    input:

    tuple path(pop_a), path(pop_b)

    path vcf

    path tbi


    output:

    path "*.weir.fst", emit: fst

    path "*.log",      emit: log


    script:

    def group_a = pop_a.name.replace("_samples.txt", "")

    def group_b = pop_b.name.replace("_samples.txt", "")

    def comparison = "${group_a}_vs_${group_b}"

    //
    // vcftools 0.1.13 -- the build on this cluster -- prints its summary
    // to stderr and writes no .log file, whatever its own documentation
    // says. Declaring "*.log" as an output without capturing it makes
    // Nextflow fail every task with MissingFileException *after*
    // vcftools has exited 0 and written a perfectly good .weir.fst, so
    // nothing is published and the whole run is marked FAILED.
    //
    // The summary is worth keeping: summarize_fst.py reads the mean and
    // weighted estimates from it and checks them against what it
    // computes from the per-SNP table. Newer builds do write the file
    // themselves, so it is only filled in when missing.
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

    """

}


process SUMMARIZE_FST {

    tag "Summarise ${params.k}-group comparisons"


    publishDir "${params.fst_outdir}/pairwise_fst",
        mode: "copy"


    input:

    path fst_files

    path log_files


    output:

    path "pairwise_fst_summary_K*.tsv"


    script:

    """

    summarize_fst.py \
        --fstdir . \
        --output pairwise_fst_summary_K${params.k}.tsv

    """

}


process PLOT_FST_HEATMAP {

    tag "FST heatmap K=${params.k}"


    publishDir "${params.fst_outdir}/plots",
        mode: "copy"


    input:

    path summary

    path labels


    output:

    path "pairwise_fst_heatmap_K*.pdf"

    path "pairwise_fst_heatmap_K*.png"

    path "pairwise_fst_matrix_K*.csv"


    script:

    """

    plot_fst_heatmap.py \
        --summary ${summary} \
        --labels ${labels} \
        --metric ${params.metric} \
        --k ${params.k} \
        --outdir .

    """

}


/*
 * K3_samples.txt -> 3
 */

def group_index(path) {

    return path.name.replaceAll(/^K(\d+)_samples\.txt$/, '$1') as Integer
}


/*
 * --pairs all          every pair
 * --pairs K1:K5        one comparison
 * --pairs K1:K5,K2:K3  several
 */

def wanted_pairs(spec) {

    if (spec == null || spec.toString().toLowerCase() == "all") {
        return null
    }

    return spec.toString()
        .split(",")
        .collect { it.trim().split("[:_-]+").findAll { it }.sort().join("_") }
        .toSet()
}


workflow {

    fam = Channel.value(
        file(params.fam, checkIfExists: true)
    )

    metadata = Channel.value(
        file(params.metadata, checkIfExists: true)
    )

    vcf = Channel.value(
        file(params.vcf, checkIfExists: true)
    )

    tbi = Channel.value(
        file("${params.vcf}.tbi", checkIfExists: true)
    )

    /*
     * Checked by hand rather than with checkIfExists, because the bare
     * "No such file or directory" is a poor answer to the most common
     * way this pipeline fails: ADMIXTURE has not been run for this K,
     * or its output is not where --admixture_dir says.
     */

    def q_path = file(params.q_file)

    if (!q_path.exists()) {

        def available = file(params.admixture_dir).exists()
            ? file("${params.admixture_dir}/*.Q").collect { it.name }.sort()
            : []

        error """
        No ADMIXTURE Q file for K=${params.k}.

          looked for : ${params.q_file}
          in         : ${params.admixture_dir}
          found there: ${available ? available.join(', ') : 'nothing, or the directory does not exist'}

        ADMIXTURE writes .Q into the working directory of the job that
        ran it, so a personal copy may be elsewhere:

          find ~ -name 'cauc_filtered.final.*.Q'

        Or point at the shared copy:

          --admixture_dir /path/to/shared-data/admixture/results
        """.stripIndent()
    }

    q_file = Channel.value(q_path)


    ASSIGN_GROUPS(
        fam,
        q_file
    )


    DESCRIBE_GROUPS(
        ASSIGN_GROUPS.out.q_values,
        metadata
    )


    /*
     * Groups too small for a meaningful FST estimate are dropped here
     * rather than left to produce a column of NaN.
     */

    usable = ASSIGN_GROUPS.out.lists
        .flatten()
        .filter { it.countLines() >= params.min_group_size }


    selected = wanted_pairs(params.pairs)


    pairs = usable
        .toSortedList { a, b -> group_index(a) <=> group_index(b) }
        .flatMap { lists ->

            def combinations = []

            for (int i = 0; i < lists.size(); i++) {
                for (int j = i + 1; j < lists.size(); j++) {
                    combinations << [ lists[i], lists[j] ]
                }
            }

            return combinations
        }
        .filter { pair ->

            if (selected == null) {
                return true
            }

            def key = pair
                .collect { it.name.replace("_samples.txt", "") }
                .sort()
                .join("_")

            return selected.contains(key)
        }


    PAIRWISE_FST(
        pairs,
        vcf,
        tbi
    )


    SUMMARIZE_FST(
        PAIRWISE_FST.out.fst.collect(),
        PAIRWISE_FST.out.log.collect()
    )


    PLOT_FST_HEATMAP(
        SUMMARIZE_FST.out,
        DESCRIBE_GROUPS.out.labels
    )

}


workflow.onComplete {

    log.info """
    ----------------------------------------------------------
    FST, K = ${params.k}
    results : ${params.fst_outdir}
    status  : ${workflow.success ? 'OK' : 'FAILED'}
    ----------------------------------------------------------
    """.stripIndent()
}
