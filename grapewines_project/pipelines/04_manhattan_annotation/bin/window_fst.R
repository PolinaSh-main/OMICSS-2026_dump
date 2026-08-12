#!/usr/bin/env Rscript

# ---------------------------------------------------------------------
#  Per-SNP FST  ->  windowed FST  ->  Manhattan plot
#
#  The windowing is done here rather than by vcftools --fst-window-size,
#  because the instructions ask for the *.weir.fst files that already
#  exist to be reused.
#
#  Each chromosome is cut into fixed windows (50 kb by default) and the
#  mean FST of the SNPs inside each window is taken. Windows carrying
#  fewer than --min-snps SNPs are dropped: a window with three SNPs can
#  reach a high mean by chance and would otherwise dominate the tail.
#  What remains is ranked, and the top --top-fraction is marked as
#  outlying.
#
#  One point on the plot is one window, not one SNP.
#
#  Usage:
#      window_fst.R --fst K1_vs_K3.weir.fst --outdir results/
# ---------------------------------------------------------------------


parse_args <- function(argv) {

    args <- list(
        fst           = NA_character_,
        outdir        = ".",
        window        = "50000",
        min_snps      = "20",
        top_fraction  = "0.01",
        comparison    = NA_character_
    )

    i <- 1

    while (i <= length(argv)) {

        key <- sub("^--", "", argv[i])

        if (!key %in% names(args)) {
            stop("Unknown argument: ", argv[i])
        }

        if (i + 1 > length(argv)) {
            stop("Missing value for ", argv[i])
        }

        args[[key]] <- argv[i + 1]

        i <- i + 2
    }

    if (is.na(args$fst)) {
        stop("--fst is required")
    }

    args$window       <- as.numeric(args$window)
    args$min_snps     <- as.integer(args$min_snps)
    args$top_fraction <- as.numeric(args$top_fraction)

    if (is.na(args$comparison)) {
        args$comparison <- sub("\\.weir\\.fst$", "", basename(args$fst))
    }

    args
}


read_fst <- function(path) {

    fst <- read.delim(
        path,
        stringsAsFactors = FALSE,
        colClasses = c("character", "integer", "character")
    )

    expected <- c("CHROM", "POS", "WEIR_AND_COCKERHAM_FST")

    if (!all(expected %in% names(fst))) {
        stop(
            "Unexpected columns in ", path, ": ",
            paste(names(fst), collapse = ", ")
        )
    }

    # vcftools writes "-nan" where FST is undefined (a site that is
    # monomorphic in both groups). as.numeric turns those into NA.
    fst$FST <- suppressWarnings(
        as.numeric(fst$WEIR_AND_COCKERHAM_FST)
    )

    fst <- fst[!is.na(fst$FST), c("CHROM", "POS", "FST")]

    if (nrow(fst) == 0) {
        stop("No usable FST values in ", path)
    }

    fst
}


# Chromosomes sort as 1, 2, ... 10, 11 -- not 1, 10, 11, 2.
chromosome_order <- function(chroms) {

    unique_chroms <- unique(chroms)

    numeric_part <- suppressWarnings(
        as.numeric(gsub("[^0-9]", "", unique_chroms))
    )

    unique_chroms[order(numeric_part, unique_chroms, na.last = TRUE)]
}


make_windows <- function(fst, window_size, min_snps) {

    # Window index, then the coordinates that index stands for.
    idx <- floor(fst$POS / window_size)

    key <- paste(fst$CHROM, idx, sep = ":")

    agg <- data.frame(
        key       = names(tapply(fst$FST, key, mean)),
        mean_fst  = as.numeric(tapply(fst$FST, key, mean)),
        max_fst   = as.numeric(tapply(fst$FST, key, max)),
        n_snps    = as.integer(tapply(fst$FST, key, length)),
        stringsAsFactors = FALSE
    )

    split_key <- do.call(rbind, strsplit(agg$key, ":", fixed = TRUE))

    windows <- data.frame(
        chrom        = split_key[, 1],
        window_start = as.numeric(split_key[, 2]) * window_size,
        stringsAsFactors = FALSE
    )

    windows$window_end <- windows$window_start + window_size
    windows$n_snps     <- agg$n_snps
    windows$mean_fst   <- agg$mean_fst
    windows$max_snp_fst <- agg$max_fst

    kept <- windows[windows$n_snps >= min_snps, ]

    if (nrow(kept) == 0) {
        stop(
            "No window has at least ", min_snps, " SNPs. ",
            "Lower --min-snps or widen --window."
        )
    }

    order_levels <- chromosome_order(kept$chrom)

    kept$chrom <- factor(kept$chrom, levels = order_levels)

    kept <- kept[order(kept$chrom, kept$window_start), ]

    rownames(kept) <- NULL

    kept
}


mark_outliers <- function(windows, top_fraction) {

    threshold <- quantile(
        windows$mean_fst,
        probs = 1 - top_fraction,
        na.rm = TRUE
    )

    windows$outlier <- windows$mean_fst >= threshold

    attr(windows, "threshold") <- as.numeric(threshold)

    windows
}


draw_manhattan <- function(windows, threshold, comparison, window_size) {

    chroms <- levels(windows$chrom)

    # Lay the chromosomes end to end along the x axis.
    widths <- tapply(
        windows$window_end,
        windows$chrom,
        max
    )

    widths[is.na(widths)] <- 0

    offsets <- c(0, cumsum(as.numeric(widths))[-length(widths)])

    names(offsets) <- chroms

    x <- offsets[as.character(windows$chrom)] + windows$window_start

    # Alternating shades so neighbouring chromosomes stay apart, with
    # the outlying windows picked out in the project's accent colour.
    base_colours <- c("#4E6E80", "#9DB4C0")

    point_colour <- base_colours[
        (match(as.character(windows$chrom), chroms) %% 2) + 1
    ]

    point_colour[windows$outlier] <- "#AE0039"

    par(mar = c(5, 5, 4, 2))

    plot(
        x, windows$mean_fst,
        pch  = 20,
        cex  = ifelse(windows$outlier, 0.9, 0.5),
        col  = point_colour,
        xaxt = "n",
        xlab = "Chromosome",
        ylab = sprintf("Mean FST per %g kb window", window_size / 1000),
        main = sprintf(
            "%s  -  %d windows, %d outliers above %.4f",
            comparison,
            nrow(windows),
            sum(windows$outlier),
            threshold
        ),
        cex.main = 1.0,
        bty = "n"
    )

    centres <- offsets + as.numeric(widths) / 2

    axis(
        1,
        at     = centres,
        labels = chroms,
        las    = 2,
        cex.axis = 0.7,
        tick   = FALSE
    )

    abline(
        h   = threshold,
        col = "#AE0039",
        lty = 2
    )
}


main <- function() {

    args <- parse_args(commandArgs(trailingOnly = TRUE))

    dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

    fst <- read_fst(args$fst)

    cat(sprintf(
        "%s: %d usable SNPs on %d chromosomes\n",
        args$comparison,
        nrow(fst),
        length(unique(fst$CHROM))
    ))

    windows <- make_windows(fst, args$window, args$min_snps)

    windows <- mark_outliers(windows, args$top_fraction)

    threshold <- attr(windows, "threshold")

    cat(sprintf(
        "%d windows with >= %d SNPs, %d outliers (top %.1f%%, FST >= %.4f)\n",
        nrow(windows),
        args$min_snps,
        sum(windows$outlier),
        args$top_fraction * 100,
        threshold
    ))


    windows_path <- file.path(
        args$outdir,
        sprintf("%s_windows.tsv", args$comparison)
    )

    write.table(
        windows,
        windows_path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )


    outliers_path <- file.path(
        args$outdir,
        sprintf("%s_outlier_windows.tsv", args$comparison)
    )

    write.table(
        windows[windows$outlier, ],
        outliers_path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )


    pdf_path <- file.path(
        args$outdir,
        sprintf("%s_manhattan.pdf", args$comparison)
    )

    pdf(pdf_path, width = 11, height = 4.5)
    draw_manhattan(windows, threshold, args$comparison, args$window)
    invisible(dev.off())

    png_path <- file.path(
        args$outdir,
        sprintf("%s_manhattan.png", args$comparison)
    )

    png(png_path, width = 11, height = 4.5, units = "in", res = 200)
    draw_manhattan(windows, threshold, args$comparison, args$window)
    invisible(dev.off())


    cat("Wrote:\n")
    cat(" ", windows_path, "\n")
    cat(" ", outliers_path, "\n")
    cat(" ", pdf_path, "\n")
    cat(" ", png_path, "\n")
}


main()
