#!/usr/bin/env python3
"""
Per-SNP FST  ->  windowed FST  ->  Manhattan plot.

The windowing happens here rather than through vcftools
--fst-window-size, because the task is to reuse the *.weir.fst files
that already exist.

Each chromosome is cut into fixed windows (50 kb by default) and the
mean FST of the SNPs inside each window is taken. Windows holding fewer
than --min-snps SNPs are dropped: a window built from three SNPs can
reach a high mean by chance and would otherwise fill the tail. What
remains is ranked and the top --top-fraction is marked as outlying.

One point on the plot is one window, not one SNP.

Output:
    <comparison>_windows.tsv
    <comparison>_outlier_windows.tsv
    <comparison>_manhattan.pdf
    <comparison>_manhattan.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


#
# Alternating shades so neighbouring chromosomes stay apart, plus the
# project's crimson for the outliers.
#

CHROMOSOME_COLOURS = ["#4E6E80", "#9DB4C0"]

OUTLIER_COLOUR = "#AE0039"


def parse_args():

    parser = argparse.ArgumentParser(
        description="Window per-SNP FST and draw a Manhattan plot"
    )

    parser.add_argument(
        "--fst",
        required=True,
        type=Path,
        help="vcftools *.weir.fst"
    )

    parser.add_argument(
        "--window",
        type=int,
        default=50_000,
        help="Window size in bp (default 50000)"
    )

    parser.add_argument(
        "--min-snps",
        type=int,
        default=20,
        help="Drop windows with fewer SNPs than this (default 20)"
    )

    parser.add_argument(
        "--top-fraction",
        type=float,
        default=0.01,
        help="Fraction of windows marked as outlying (default 0.01)"
    )

    parser.add_argument(
        "--comparison",
        default=None,
        help="Name used in the output files (default: from --fst)"
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path
    )

    return parser.parse_args()


def chromosome_sort_key(value):
    """
    1, 2, ... 10 rather than 1, 10, 2. Anything without digits sorts
    last, alphabetically.
    """

    digits = "".join(c for c in str(value) if c.isdigit())

    return (0, int(digits), "") if digits else (1, 0, str(value))


def read_fst(path: Path) -> pd.DataFrame:

    fst = pd.read_csv(
        path,
        sep="\t",
        dtype={"CHROM": str, "POS": "int64"},
        #
        # vcftools writes "-nan" where FST is undefined, which happens
        # at sites that are monomorphic in both groups.
        #
        na_values=["-nan", "nan", "NaN", "inf", "-inf"]
    )

    expected = {"CHROM", "POS", "WEIR_AND_COCKERHAM_FST"}

    missing = expected - set(fst.columns)

    if missing:
        raise ValueError(
            f"{path} is missing columns {sorted(missing)}; "
            f"found {list(fst.columns)}"
        )

    fst = fst.rename(columns={"WEIR_AND_COCKERHAM_FST": "fst"})

    fst = fst.dropna(subset=["fst"])

    if fst.empty:
        raise ValueError(f"No usable FST values in {path}")

    return fst[["CHROM", "POS", "fst"]]


def make_windows(
    fst: pd.DataFrame,
    window_size: int,
    min_snps: int
) -> pd.DataFrame:

    fst = fst.assign(
        window_start=(fst["POS"] // window_size) * window_size
    )

    grouped = fst.groupby(
        ["CHROM", "window_start"],
        sort=False
    )["fst"]

    windows = grouped.agg(
        n_snps="size",
        mean_fst="mean",
        max_snp_fst="max"
    ).reset_index()

    windows = windows.rename(columns={"CHROM": "chrom"})

    windows["window_end"] = windows["window_start"] + window_size


    kept = windows[windows["n_snps"] >= min_snps].copy()

    if kept.empty:
        raise ValueError(
            f"No window holds at least {min_snps} SNPs. "
            f"Lower --min-snps or widen --window."
        )

    kept = (
        kept
        .assign(_chrom_key=kept["chrom"].map(chromosome_sort_key))
        .sort_values(["_chrom_key", "window_start"])
        .drop(columns="_chrom_key")
        .reset_index(drop=True)
    )

    return kept[
        [
            "chrom",
            "window_start",
            "window_end",
            "n_snps",
            "mean_fst",
            "max_snp_fst",
        ]
    ]


def mark_outliers(
    windows: pd.DataFrame,
    top_fraction: float
) -> tuple[pd.DataFrame, float]:

    threshold = float(
        windows["mean_fst"].quantile(1.0 - top_fraction)
    )

    windows = windows.assign(
        outlier=windows["mean_fst"] >= threshold
    )

    return windows, threshold


def draw_manhattan(
    windows: pd.DataFrame,
    threshold: float,
    comparison: str,
    window_size: int
):

    chromosomes = sorted(
        windows["chrom"].unique(),
        key=chromosome_sort_key
    )

    #
    # Lay the chromosomes end to end along the x axis.
    #

    widths = (
        windows.groupby("chrom")["window_end"]
        .max()
        .reindex(chromosomes)
        .fillna(0)
    )

    offsets = widths.cumsum().shift(fill_value=0)

    x = (
        windows["chrom"].map(offsets)
        + windows["window_start"]
    ).to_numpy(dtype=float)

    shade = np.array(
        [
            CHROMOSOME_COLOURS[chromosomes.index(c) % 2]
            for c in windows["chrom"]
        ]
    )

    outlier = windows["outlier"].to_numpy()

    figure, axes = plt.subplots(figsize=(11, 4.5))

    axes.scatter(
        x[~outlier],
        windows.loc[~outlier, "mean_fst"],
        s=5,
        c=shade[~outlier],
        linewidths=0,
        rasterized=True
    )

    axes.scatter(
        x[outlier],
        windows.loc[outlier, "mean_fst"],
        s=14,
        c=OUTLIER_COLOUR,
        linewidths=0,
        rasterized=True
    )

    axes.axhline(
        threshold,
        color=OUTLIER_COLOUR,
        linestyle="--",
        linewidth=0.9
    )

    centres = offsets + widths / 2.0

    axes.set_xticks(centres.to_numpy())
    axes.set_xticklabels(chromosomes, fontsize=8)

    axes.set_xlabel("Chromosome")

    axes.set_ylabel(
        f"Mean FST per {window_size / 1000:g} kb window"
    )

    axes.set_title(
        f"{comparison}  -  {len(windows)} windows, "
        f"{int(outlier.sum())} outliers above {threshold:.4f}",
        fontsize=11
    )

    axes.set_xlim(-0.01 * x.max(), 1.01 * x.max())

    axes.spines["top"].set_visible(False)
    axes.spines["right"].set_visible(False)

    axes.tick_params(length=0)

    figure.tight_layout()

    return figure


def main():

    args = parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    comparison = args.comparison or args.fst.name.replace(
        ".weir.fst", ""
    )


    fst = read_fst(args.fst)

    print(
        f"{comparison}: {len(fst)} usable SNPs on "
        f"{fst['CHROM'].nunique()} chromosomes"
    )


    windows = make_windows(fst, args.window, args.min_snps)

    windows, threshold = mark_outliers(windows, args.top_fraction)

    print(
        f"{len(windows)} windows with >= {args.min_snps} SNPs, "
        f"{int(windows['outlier'].sum())} outliers "
        f"(top {args.top_fraction:.1%}, FST >= {threshold:.4f})"
    )


    windows_path = args.outdir / f"{comparison}_windows.tsv"

    windows.to_csv(
        windows_path,
        sep="\t",
        index=False,
        float_format="%.6f"
    )


    outliers_path = args.outdir / f"{comparison}_outlier_windows.tsv"

    windows[windows["outlier"]].to_csv(
        outliers_path,
        sep="\t",
        index=False,
        float_format="%.6f"
    )


    figure = draw_manhattan(
        windows,
        threshold,
        comparison,
        args.window
    )

    pdf_path = args.outdir / f"{comparison}_manhattan.pdf"
    png_path = args.outdir / f"{comparison}_manhattan.png"

    figure.savefig(pdf_path)
    figure.savefig(png_path, dpi=200)

    plt.close(figure)


    print("Wrote:")
    print(" ", windows_path)
    print(" ", outliers_path)
    print(" ", pdf_path)
    print(" ", png_path)


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
