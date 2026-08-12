#!/usr/bin/env Rscript

# Bootstrap the SNPhylo maximum-likelihood tree.
#
# This replaces SNPhylo's determine_bs_tree.R for one reason: that script
# calls
#
#     bootstrap.pml(fit, bs = n, optNni = TRUE, multicore = TRUE)
#
# without mc.cores, so phangorn falls back to parallel::detectCores().
# On these compute nodes that reports 64 regardless of what SLURM
# allocated, so a job asking for 5 CPUs forks 64 workers, each holding
# its own copy of a 413 x 12848 alignment. The group's earlier attempt
# was killed that way. Here the worker count is passed in and matches
# the allocation.
#
# Everything else -- the number of replicates, the NNI optimisation, the
# seed -- is left as SNPhylo had it.

suppressPackageStartupMessages({
    library(getopt)
    library(phangorn)
})

spec <- matrix(c(
    "tree",       "i", 1, "character",
    "phylip",     "p", 1, "character",
    "out",        "o", 1, "character",
    "replicates", "n", 1, "integer",
    "cores",      "c", 1, "integer",
    "help",       "h", 0, "logical"
), ncol = 4, byrow = TRUE)

opt <- getopt(spec)

if (!is.null(opt$help) ||
    is.null(opt$tree) || is.null(opt$phylip) || is.null(opt$out)) {

    cat("Usage: bootstrap_tree.R --tree ml.tree --phylip aln.phylip",
        "--out bs.tree [--replicates 100] [--cores 8]\n")

    quit(save = "no", status = if (is.null(opt$help)) 1 else 0)
}

replicates <- if (is.null(opt$replicates)) 100L else opt$replicates
cores      <- if (is.null(opt$cores)) 1L else opt$cores

stopifnot(file.exists(opt$tree), file.exists(opt$phylip))

alignment <- read.phyDat(opt$phylip, format = "interleaved", type = "DNA")
ml_tree   <- read.tree(opt$tree)

cat(sprintf("%d tips, %d sites, %d replicates on %d cores\n",
            length(ml_tree$tip.label), sum(attr(alignment, "weight")),
            replicates, cores))

fit <- pml(ml_tree, alignment)

set.seed(1)

replicate_trees <- bootstrap.pml(
    fit,
    bs        = replicates,
    optNni    = TRUE,
    multicore = cores > 1,
    mc.cores  = cores
)

# plotBS() is what maps the replicate topologies onto the reference tree
# as node labels. The device is a throwaway: the figures come from
# plot_tree.py, this call is only here for its return value.
grDevices::pdf(file = NULL)
options(warn = -1)
bs_tree <- plotBS(fit$tree, replicate_trees, type = "none")
options(warn = 0)
invisible(grDevices::dev.off())

# phangorn 2.12 returns these as proportions in [0, 1]. Every convention
# for reading a tree -- and the "support >= 70" the task asks for -- is
# in percent, so rescale before writing rather than leaving a file whose
# numbers mean something different from what they look like.

support <- suppressWarnings(as.numeric(bs_tree$node.label))

finite <- support[is.finite(support)]

if (length(finite) > 0 && max(finite) <= 1) {

    support <- support * 100

    cat("node labels were proportions; rescaled to percent\n")
}

bs_tree$node.label <- ifelse(
    is.finite(support),
    format(round(support), trim = TRUE, scientific = FALSE),
    ""
)

write.tree(bs_tree, file = opt$out)

cat(sprintf("wrote %s\n", opt$out))

cat(sprintf("support >= 70: %d of %d internal nodes\n",
            sum(finite * (if (max(finite) <= 1) 100 else 1) >= 70),
            length(finite)))
