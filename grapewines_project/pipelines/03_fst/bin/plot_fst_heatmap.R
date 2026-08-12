#!/usr/bin/env Rscript

# ---------------------------------------------------------------------
#  Pairwise FST heatmap.
#
#  Reads the summary table produced by summarize_fst.py, folds it into a
#  square matrix and draws it. Lighter cells are genetically similar
#  groups, darker cells are differentiated ones.
#
#  Written against base R only. The cluster has no guaranteed ggplot2,
#  and a missing package at the end of a 21-job FST run is an expensive
#  way to find out.
#
#  Usage:
#      plot_fst_heatmap.R --summary <tsv> --outdir <dir> [options]
# ---------------------------------------------------------------------


parse_args <- function(argv) {

    defaults <- list(
        summary = NA_character_,
        labels  = NA_character_,
        outdir  = ".",
        metric  = "mean_fst",
        k       = NA_character_
    )

    i <- 1

    while (i <= length(argv)) {

        key <- sub("^--", "", argv[i])

        if (!key %in% names(defaults)) {
            stop("Unknown argument: ", argv[i])
        }

        if (i + 1 > length(argv)) {
            stop("Missing value for ", argv[i])
        }

        defaults[[key]] <- argv[i + 1]

        i <- i + 2
    }

    if (is.na(defaults$summary)) {
        stop("--summary is required")
    }

    defaults
}


natural_group_order <- function(groups) {

    # K10 must sort after K9, not after K1.

    groups[order(as.integer(sub("^K", "", groups)))]
}


build_matrix <- function(summary, metric) {

    if (!metric %in% names(summary)) {
        stop(
            "Column '", metric, "' not in the summary table. ",
            "Available: ", paste(names(summary), collapse = ", ")
        )
    }

    groups <- natural_group_order(
        unique(c(summary$group_a, summary$group_b))
    )

    m <- matrix(
        NA_real_,
        nrow = length(groups),
        ncol = length(groups),
        dimnames = list(groups, groups)
    )

    for (row in seq_len(nrow(summary))) {

        a <- summary$group_a[row]
        b <- summary$group_b[row]
        v <- summary[[metric]][row]

        m[a, b] <- v
        m[b, a] <- v
    }

    # A group has no differentiation from itself.
    diag(m) <- 0

    m
}


# Lines of margin needed to fit the axis labels, at cex.axis = 0.8.
label_margin <- function(labels) {

    min(20, max(5, max(nchar(labels)) * 0.45))
}


draw_heatmap <- function(m, labels, metric, main_title) {

    n <- nrow(m)

    palette <- colorRampPalette(
        c("#FBF7F0", "#BCDCE8", "#4E9DBC", "#00699A", "#0B3A4F")
    )(256)

    off_diagonal <- m[upper.tri(m)]

    if (all(is.na(off_diagonal))) {
        stop("Every off-diagonal FST value is missing; nothing to plot")
    }

    limits <- c(0, max(off_diagonal, na.rm = TRUE))

    margin <- label_margin(labels)

    par(mar = c(margin, margin, 4, 7))

    image(
        x    = seq_len(n),
        y    = seq_len(n),
        z    = t(m[n:1, , drop = FALSE]),
        col  = palette,
        zlim = limits,
        axes = FALSE,
        xlab = "",
        ylab = ""
    )

    axis(
        1,
        at     = seq_len(n),
        labels = labels,
        las    = 2,
        cex.axis = 0.8,
        tick   = FALSE
    )

    axis(
        2,
        at     = seq_len(n),
        labels = rev(labels),
        las    = 2,
        cex.axis = 0.8,
        tick   = FALSE
    )

    # Cell values. White text on the dark end so it stays readable.
    for (i in seq_len(n)) {
        for (j in seq_len(n)) {

            value <- m[n - j + 1, i]

            if (is.na(value)) next

            shade <- (value - limits[1]) / diff(limits)

            text(
                i, j,
                labels = formatC(value, format = "f", digits = 4),
                cex    = 0.7,
                col    = if (shade > 0.55) "white" else "grey15"
            )
        }
    }

    box(col = "grey80")

    title(main = main_title, cex.main = 1.1)

    # Colour key in the right-hand margin.
    key_at <- seq(limits[1], limits[2], length.out = 5)

    legend(
        "topright",
        inset  = c(-0.09, 0),
        legend = formatC(rev(key_at), format = "f", digits = 3),
        fill   = rev(palette[round(seq(1, 256, length.out = 5))]),
        border = NA,
        bty    = "n",
        xpd    = TRUE,
        cex    = 0.7,
        title  = metric
    )
}


main <- function() {

    args <- parse_args(commandArgs(trailingOnly = TRUE))

    summary <- read.delim(
        args$summary,
        stringsAsFactors = FALSE
    )

    m <- build_matrix(summary, args$metric)

    groups <- rownames(m)

    k <- if (is.na(args$k)) length(groups) else args$k


    # Axis labels: interpreted names when describe_groups.py provided
    # them, bare K<i> otherwise.

    labels <- groups

    if (!is.na(args$labels) && file.exists(args$labels)) {

        lookup <- read.delim(args$labels, stringsAsFactors = FALSE)

        matched <- match(groups, lookup$group)

        labels <- ifelse(
            is.na(matched),
            groups,
            lookup$label[matched]
        )
    }


    dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

    matrix_path <- file.path(
        args$outdir,
        sprintf("pairwise_fst_matrix_K%s.csv", k)
    )

    write.csv(
        m,
        matrix_path,
        quote = FALSE
    )


    heading <- sprintf(
        "Pairwise FST between ADMIXTURE groups (K = %s)",
        k
    )

    # Figure size = plot area + whatever the labels need. One line of
    # margin is 0.2 inch at the default point size.

    side <- 2 + 0.6 * length(groups) + 0.2 * label_margin(labels)

    pdf_path <- file.path(
        args$outdir,
        sprintf("pairwise_fst_heatmap_K%s.pdf", k)
    )

    pdf(pdf_path, width = side, height = side)
    draw_heatmap(m, labels, args$metric, heading)
    invisible(dev.off())

    png_path <- file.path(
        args$outdir,
        sprintf("pairwise_fst_heatmap_K%s.png", k)
    )

    png(png_path, width = side, height = side, units = "in", res = 200)
    draw_heatmap(m, labels, args$metric, heading)
    invisible(dev.off())


    cat("Wrote:\n")
    cat(" ", matrix_path, "\n")
    cat(" ", pdf_path, "\n")
    cat(" ", png_path, "\n")

    cat("\nLowest  FST: ")
    off <- m
    diag(off) <- NA
    lo <- which(off == min(off, na.rm = TRUE), arr.ind = TRUE)[1, ]
    cat(rownames(m)[lo[1]], "vs", colnames(m)[lo[2]],
        "=", formatC(min(off, na.rm = TRUE), format = "f", digits = 4), "\n")

    cat("Highest FST: ")
    hi <- which(off == max(off, na.rm = TRUE), arr.ind = TRUE)[1, ]
    cat(rownames(m)[hi[1]], "vs", colnames(m)[hi[2]],
        "=", formatC(max(off, na.rm = TRUE), format = "f", digits = 4), "\n")
}


main()
