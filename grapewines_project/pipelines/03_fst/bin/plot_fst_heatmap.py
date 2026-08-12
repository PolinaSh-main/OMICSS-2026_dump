#!/usr/bin/env python3
"""
Pairwise FST heatmap.

Reads the summary table from summarize_fst.py, folds it into a square
matrix and draws it. Light cells are groups that are genetically
similar, dark cells are differentiated ones.

Input:
    pairwise_fst_summary_K<K>.tsv
    group_labels_K<K>.tsv          optional, from describe_groups.py

Output:
    pairwise_fst_heatmap_K<K>.pdf
    pairwise_fst_heatmap_K<K>.png
    pairwise_fst_matrix_K<K>.csv
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
from matplotlib.colors import LinearSegmentedColormap


#
# Sequential ramp built around the project's K1 blue, so the FST
# figures sit next to the ADMIXTURE barplots without clashing.
#

FST_COLOURS = [
    "#FBF7F0",
    "#BCDCE8",
    "#4E9DBC",
    "#00699A",
    "#0B3A4F",
]


def parse_args():

    parser = argparse.ArgumentParser(
        description="Heatmap of pairwise FST between ADMIXTURE groups"
    )

    parser.add_argument(
        "--summary",
        required=True,
        type=Path,
        help="pairwise_fst_summary_K<K>.tsv"
    )

    parser.add_argument(
        "--labels",
        type=Path,
        default=None,
        help="group_labels_K<K>.tsv; falls back to bare K<i>"
    )

    parser.add_argument(
        "--metric",
        default="mean_fst",
        help="Column of the summary table to plot (default mean_fst)"
    )

    parser.add_argument(
        "--k",
        default=None,
        help="K, used in the output file names "
             "(default: number of groups found)"
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path
    )

    return parser.parse_args()


def group_number(group: str) -> int:
    """
    K10 must sort after K9, not after K1.
    """

    return int(str(group).lstrip("K"))


def build_matrix(summary: pd.DataFrame, metric: str) -> pd.DataFrame:

    if metric not in summary.columns:
        raise ValueError(
            f"Column '{metric}' is not in the summary table. "
            f"Available: {list(summary.columns)}"
        )

    groups = sorted(
        set(summary["group_a"]) | set(summary["group_b"]),
        key=group_number
    )

    position = {group: i for i, group in enumerate(groups)}

    #
    # Filled as a plain array rather than through DataFrame.loc: recent
    # pandas hands out a read-only view from .values, so writing to the
    # diagonal in place is not portable across versions.
    #

    values = np.full((len(groups), len(groups)), np.nan)

    for row in summary.itertuples():

        i = position[row.group_a]
        j = position[row.group_b]

        value = float(getattr(row, metric))

        values[i, j] = value
        values[j, i] = value

    # A group is not differentiated from itself.
    np.fill_diagonal(values, 0.0)

    return pd.DataFrame(values, index=groups, columns=groups)


def read_labels(path: Path | None, groups: list[str]) -> list[str]:

    if path is None or not path.exists():
        return list(groups)

    lookup = pd.read_csv(path, sep="\t")

    if not {"group", "label"} <= set(lookup.columns):
        print(
            f"WARNING: {path.name} has no group/label columns; "
            f"using bare group names",
            file=sys.stderr
        )
        return list(groups)

    mapping = dict(zip(lookup["group"], lookup["label"]))

    return [mapping.get(group, group) for group in groups]


def draw(matrix: pd.DataFrame, labels: list[str], metric: str, title: str):

    n = len(matrix)

    values = matrix.to_numpy(dtype=float).copy()

    diagonal = np.eye(n, dtype=bool)

    off_diagonal = values[~diagonal]

    if np.all(np.isnan(off_diagonal)):
        raise ValueError(
            "Every off-diagonal FST value is missing; nothing to plot"
        )

    #
    # The diagonal is zero by construction and carries no information,
    # but leaving it in would stretch the colour scale from 0 and push
    # every real value into the same dark shade. Mask it, and scale the
    # ramp across the range that is actually there. Every cell is
    # labelled with its value, and the colour bar is annotated, so the
    # relative scale cannot be read as an absolute one.
    #

    values[diagonal] = np.nan

    vmin = float(np.nanmin(off_diagonal))

    vmax = float(np.nanmax(off_diagonal))

    if vmax - vmin < 1e-9:
        vmin, vmax = vmin - 1e-6, vmax + 1e-6

    colours = LinearSegmentedColormap.from_list("fst", FST_COLOURS)

    colours.set_bad("#EFEBE5")


    # Room for the labels, which carry the interpreted group names and
    # are therefore long.
    longest = max(len(label) for label in labels)

    side = 2.4 + 0.62 * n + 0.075 * longest

    figure, axes = plt.subplots(figsize=(side, side))

    image = axes.imshow(
        np.ma.masked_invalid(values),
        cmap=colours,
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest"
    )


    axes.set_xticks(range(n))
    axes.set_yticks(range(n))

    axes.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    axes.set_yticklabels(labels, fontsize=9)

    axes.set_xticks(np.arange(n + 1) - 0.5, minor=True)
    axes.set_yticks(np.arange(n + 1) - 0.5, minor=True)

    axes.grid(which="minor", color="white", linewidth=1.5)

    axes.tick_params(which="both", length=0)

    for spine in axes.spines.values():
        spine.set_visible(False)


    #
    # Print the value in every cell. The text has to flip to white on
    # the dark end of the ramp to stay readable.
    #

    for i in range(n):
        for j in range(n):

            value = values[i, j]

            if np.isnan(value):
                continue

            shade = (value - vmin) / (vmax - vmin)

            axes.text(
                j, i,
                f"{value:.4f}",
                ha="center",
                va="center",
                fontsize=8,
                color="white" if shade > 0.55 else "#222222"
            )


    axes.set_title(title, fontsize=11, pad=14)

    bar = figure.colorbar(image, ax=axes, shrink=0.72, pad=0.02)

    bar.set_label(metric, fontsize=9)

    bar.ax.tick_params(labelsize=8)

    bar.outline.set_visible(False)

    figure.tight_layout()

    return figure


def report_extremes(matrix: pd.DataFrame):

    values = matrix.to_numpy(dtype=float).copy()

    np.fill_diagonal(values, np.nan)

    if np.all(np.isnan(values)):
        return

    lowest = np.unravel_index(np.nanargmin(values), values.shape)
    highest = np.unravel_index(np.nanargmax(values), values.shape)

    groups = list(matrix.index)

    print(
        f"Lowest  FST: {groups[lowest[0]]} vs {groups[lowest[1]]} "
        f"= {values[lowest]:.4f}"
    )

    print(
        f"Highest FST: {groups[highest[0]]} vs {groups[highest[1]]} "
        f"= {values[highest]:.4f}"
    )


def main():

    args = parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)


    summary = pd.read_csv(args.summary, sep="\t")

    matrix = build_matrix(summary, args.metric)

    groups = list(matrix.index)

    k = args.k if args.k is not None else len(groups)

    labels = read_labels(args.labels, groups)


    matrix_path = args.outdir / f"pairwise_fst_matrix_K{k}.csv"

    matrix.to_csv(matrix_path, float_format="%.6f")


    title = f"Pairwise FST between ADMIXTURE groups (K = {k})"

    figure = draw(matrix, labels, args.metric, title)

    pdf_path = args.outdir / f"pairwise_fst_heatmap_K{k}.pdf"
    png_path = args.outdir / f"pairwise_fst_heatmap_K{k}.png"

    figure.savefig(pdf_path)
    figure.savefig(png_path, dpi=200)

    plt.close(figure)


    print("Wrote:")
    print(" ", matrix_path)
    print(" ", pdf_path)
    print(" ", png_path)
    print()

    report_extremes(matrix)


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
