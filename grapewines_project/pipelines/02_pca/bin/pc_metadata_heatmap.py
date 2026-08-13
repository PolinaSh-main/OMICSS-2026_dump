#!/usr/bin/env python3
"""
How much of each principal component is explained by each metadata
variable.

For every PC x metadata pair the samples are split into the categories of
that variable and a one-way ANOVA is fitted on the PC scores:

    PC_score ~ category

The cell value is the adjusted R2 of that fit -- the share of the
variation along that PC that the variable accounts for, penalised for the
number of categories. Stars mark the BH-adjusted ANOVA p-value.

This is a port of the group's heatmap.r. R lives on this cluster without
ggplot2 (and without a working png device), and the rest of the project
is Python already, so the figure is drawn with matplotlib like the
ADMIXTURE and FST ones.

Input:
    <prefix>.eigenvec              PLINK PCA, no header: FID IID PC1 PC2 ...
    cauc_grape_metadata.csv        sample IDs in "Column 1"

Output:
    pc_metadata_association.png / .pdf
    pc_metadata_association.csv    every pair, including the untested ones
"""

from __future__ import annotations

import argparse
import sys
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import f as f_distribution


#
# The same ramp as the FST heatmap, so the two matrices in the deck read
# as one pair of figures. Low = pale, high = the project's K1 blue and
# darker.
#

HEATMAP_COLOURS = [
    "#FBF7F0",
    "#BCDCE8",
    "#4E9DBC",
    "#00699A",
    "#0B3A4F",
]


#
# Column of the metadata CSV holding the sample ID, and the variables
# worth testing. "Variety ID" is deliberately absent: it is unique per
# sample, so every category would have one member and nothing would
# survive the minimum-group-size filter.
#

ID_COLUMN = "Column 1"

DEFAULT_VARIABLES = [
    "Genetic background",
    "Utilization",
    "Berry skin color",
    "Flower phenotype",
    "Geographic origin by region",
    "Bunch density",
    "Geographic origin by country",
    "Muscat taste",
]


#
# Values that mean "not recorded" rather than a category of their own.
# Compared case-insensitively after trimming. Empty strings and the CSV's
# own NA markers are dropped before this.
#

MISSING_LABELS = {"UNKNOWN", "NA", "N/A", "NONE RECORDED", "-", "?"}


def parse_args():

    parser = argparse.ArgumentParser(
        description="Heatmap of PC / metadata associations (adjusted R2)"
    )

    parser.add_argument(
        "--eigenvec",
        required=True,
        type=Path,
        help="PLINK .eigenvec (no header)"
    )

    parser.add_argument(
        "--metadata",
        required=True,
        type=Path,
        help="cauc_grape_metadata.csv"
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path
    )

    parser.add_argument(
        "--n-pcs",
        type=int,
        default=5,
        help="How many leading PCs to test (default 5). The later PCs "
             "carry little variance and only made the figure wider."
    )

    parser.add_argument(
        "--min-group-size",
        type=int,
        default=3,
        help="Categories with fewer samples than this are dropped "
             "before the ANOVA (default 3)"
    )

    parser.add_argument(
        "--variables",
        nargs="+",
        default=None,
        help="Metadata columns to test; default is the eight the group "
             "settled on"
    )

    parser.add_argument(
        "--prefix",
        default="pc_metadata_association",
        help="Base name of the output files"
    )

    parser.add_argument(
        "--font-scale",
        type=float,
        default=1.0,
        help="Multiplies every font size and the cell dimensions. The "
             "default is already sized for a slide; use ~0.6 for a "
             "figure that goes into a document."
    )

    return parser.parse_args()


# ---------------------------------------------------------------------
#  Reading and matching
# ---------------------------------------------------------------------


def read_eigenvec(path: Path) -> pd.DataFrame:
    """
    PLINK writes FID, IID and then one column per PC, with no header.
    """

    table = pd.read_csv(path, sep=r"\s+", header=None)

    if table.shape[1] < 3:
        raise ValueError(
            f"{path.name} has {table.shape[1]} columns; expected "
            f"FID, IID and at least one PC"
        )

    n_pcs = table.shape[1] - 2

    table.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, n_pcs + 1)]

    return table


def clean_id(values: pd.Series) -> pd.Series:

    return values.astype(str).str.strip().str.upper()


def join_metadata(
    eigenvec: pd.DataFrame,
    metadata: pd.DataFrame,
    variables: list[str],
) -> pd.DataFrame:
    """
    Attach the metadata columns to the PCA table, matching on sample ID.

    The join is left from the PCA side: a sample without metadata keeps
    its coordinates and gets NaN categories, which the per-test filter
    then drops. That way one unmatched sample cannot silently shrink
    every other test.
    """

    if ID_COLUMN not in metadata.columns:
        raise ValueError(
            f"The metadata has no '{ID_COLUMN}' column. "
            f"Found: {list(metadata.columns)}"
        )

    missing = [v for v in variables if v not in metadata.columns]

    if missing:
        raise ValueError(
            "These metadata columns were not found: "
            + ", ".join(missing)
            + f"\nAvailable: {list(metadata.columns)}"
        )

    keyed = metadata.copy()

    keyed["_sample_id"] = clean_id(keyed[ID_COLUMN])

    duplicated = keyed["_sample_id"].duplicated().sum()

    if duplicated:
        print(
            f"WARNING: {duplicated} duplicate sample ID(s) in the "
            f"metadata; keeping the first row of each",
            file=sys.stderr
        )

        keyed = keyed.drop_duplicates("_sample_id", keep="first")

    joined = eigenvec.copy()

    joined["_sample_id"] = clean_id(joined["IID"])

    joined = joined.merge(
        keyed[["_sample_id"] + variables],
        on="_sample_id",
        how="left"
    )

    matched = joined[variables].notna().any(axis=1).sum()

    print(f"Matched {matched} of {len(joined)} PCA samples to metadata")

    if matched < len(joined):

        unmatched = joined.loc[
            ~joined[variables].notna().any(axis=1),
            "IID"
        ]

        print("Unmatched:", ", ".join(map(str, unmatched)), file=sys.stderr)

    return joined


# ---------------------------------------------------------------------
#  The test
# ---------------------------------------------------------------------


def test_association(
    scores: pd.Series,
    groups: pd.Series,
    min_group_size: int,
) -> dict:
    """
    One-way ANOVA of the PC scores across the categories of one variable.

    Returns adjusted R2, the overall F-test p-value, and how much data
    the test actually saw. Everything is NaN when there is not enough
    left to compare.
    """

    labels = groups.astype("string").str.strip()

    usable = (
        scores.notna()
        & labels.notna()
        & (labels != "")
        & ~labels.str.upper().isin(MISSING_LABELS)
    )

    scores = scores[usable].astype(float)
    labels = labels[usable]

    # Categories too small to say anything about.
    counts = labels.value_counts()

    keep = counts[counts >= min_group_size].index

    scores = scores[labels.isin(keep)]
    labels = labels[labels.isin(keep)]

    n_samples = len(scores)
    n_groups = labels.nunique()

    empty = {
        "adjusted_r2": np.nan,
        "association_strength": np.nan,
        "p_value": np.nan,
        "n_samples": n_samples,
        "n_groups": n_groups,
    }

    # A comparison needs two categories, and the residual degrees of
    # freedom (n - k) must be positive or the fit is saturated.
    if n_groups < 2 or n_samples - n_groups < 1:
        return empty

    grand_mean = scores.mean()

    group_means = scores.groupby(labels, observed=True).agg(["mean", "count"])

    ss_between = float(
        (group_means["count"] * (group_means["mean"] - grand_mean) ** 2).sum()
    )

    ss_total = float(((scores - grand_mean) ** 2).sum())

    ss_within = ss_total - ss_between

    if ss_total <= 0:
        return empty

    df_between = n_groups - 1
    df_within = n_samples - n_groups

    mean_sq_within = ss_within / df_within

    r_squared = ss_between / ss_total

    #
    # Adjusted R2 can go negative when a variable explains less than
    # chance would given its number of categories. That is a real result
    # and is kept in the CSV, but the colour scale starts at zero, so the
    # plotted value is clipped -- exactly what the R version did.
    #

    adjusted_r2 = 1.0 - (1.0 - r_squared) * (n_samples - 1) / df_within

    if mean_sq_within <= 0:
        # Every sample in a category shares a score exactly. No spread
        # left to test against.
        p_value = np.nan
    else:
        f_statistic = (ss_between / df_between) / mean_sq_within
        p_value = float(f_distribution.sf(f_statistic, df_between, df_within))

    return {
        "adjusted_r2": adjusted_r2,
        "association_strength": max(0.0, adjusted_r2),
        "p_value": p_value,
        "n_samples": n_samples,
        "n_groups": n_groups,
    }


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    """
    BH-adjusted p-values, matching R's p.adjust(method = "BH").

    NaNs stay NaN and, as in R, do not count towards the number of
    tests.
    """

    adjusted = pd.Series(np.nan, index=p_values.index, dtype=float)

    present = p_values.dropna().sort_values()

    n = len(present)

    if n == 0:
        return adjusted

    ranks = np.arange(1, n + 1)

    raw = present.to_numpy() * n / ranks

    # Enforce monotonicity from the largest p-value downwards, then cap
    # at 1.
    stepped = np.minimum.accumulate(raw[::-1])[::-1]

    adjusted.loc[present.index] = np.minimum(stepped, 1.0)

    return adjusted


def run_tests(
    data: pd.DataFrame,
    pcs: list[str],
    variables: list[str],
    min_group_size: int,
) -> pd.DataFrame:

    rows = []

    for variable in variables:
        for pc in pcs:

            result = test_association(
                data[pc],
                data[variable],
                min_group_size
            )

            rows.append({"metadata": variable, "pc": pc, **result})

    results = pd.DataFrame(rows)

    results["fdr"] = benjamini_hochberg(results["p_value"])

    results["significance"] = ""

    for threshold, stars in ((0.05, "*"), (0.01, "**"), (0.001, "***")):
        results.loc[results["fdr"] < threshold, "significance"] = stars

    results["cell_label"] = np.where(
        results["association_strength"].isna(),
        "",
        results["association_strength"].map("{:.2f}".format)
        + results["significance"]
    )

    return results


# ---------------------------------------------------------------------
#  The figure
# ---------------------------------------------------------------------


def wrap(label: str, width: int) -> str:

    return "\n".join(textwrap.wrap(label, width=width)) or label


def draw(
    results: pd.DataFrame,
    pcs: list[str],
    variables: list[str],
    min_group_size: int,
    font_scale: float,
):
    """
    PCs down the rows, metadata across the columns.

    The transpose of the R version: with only five PCs left, PCs as
    columns gave a tall narrow block, and a slide is wide.
    """

    strength = (
        results
        .pivot(index="pc", columns="metadata", values="association_strength")
        .reindex(index=pcs, columns=variables)
    )

    cell_labels = (
        results
        .pivot(index="pc", columns="metadata", values="cell_label")
        .reindex(index=pcs, columns=variables)
        .fillna("")
    )

    values = strength.to_numpy(dtype=float)

    if np.all(np.isnan(values)):
        raise ValueError(
            "Every cell is missing -- no PC/metadata pair had two "
            "categories with enough samples to compare"
        )

    n_rows, n_columns = values.shape

    vmax = float(np.nanmax(values))

    if vmax <= 0:
        # Nothing rose above zero; keep the ramp from collapsing.
        vmax = 1e-6


    # Sizes are in points and scale together, so --font-scale changes
    # how big the figure is on screen without rearranging it.
    tick_size = 23 * font_scale
    cell_size = 22 * font_scale
    title_size = 30 * font_scale
    caption_size = 15 * font_scale
    legend_size = 19 * font_scale

    wrap_width = 12

    column_labels = [wrap(v, wrap_width) for v in variables]

    label_lines = max(label.count("\n") + 1 for label in column_labels)


    #
    # Figure size is derived rather than fixed: the cells keep the same
    # physical size whatever --font-scale is, and the margins grow with
    # the labels that have to fit in them.
    #

    cell_width = 2.05 * font_scale
    cell_height = 1.45 * font_scale

    left_margin = 1.7 * font_scale
    colourbar_margin = 2.3 * font_scale
    title_margin = 1.15 * font_scale

    bottom_margin = (
        label_lines * tick_size * 1.35 / 72.0
        + 0.9 * font_scale
    )

    width = n_columns * cell_width + left_margin + colourbar_margin
    height = n_rows * cell_height + title_margin + bottom_margin

    figure, axes = plt.subplots(figsize=(width, height), layout="constrained")

    colours = LinearSegmentedColormap.from_list("pc_metadata", HEATMAP_COLOURS)

    colours.set_bad("#EFEBE5")

    image = axes.imshow(
        np.ma.masked_invalid(values),
        cmap=colours,
        vmin=0.0,
        vmax=vmax,
        aspect="auto",
        interpolation="nearest"
    )


    axes.set_xticks(range(n_columns))
    axes.set_yticks(range(n_rows))

    axes.set_xticklabels(column_labels, fontsize=tick_size)
    axes.set_yticklabels(pcs, fontsize=tick_size, fontweight="bold")

    # White gutters between the cells, drawn as a minor grid so that the
    # cells themselves stay flat colour.
    axes.set_xticks(np.arange(n_columns + 1) - 0.5, minor=True)
    axes.set_yticks(np.arange(n_rows + 1) - 0.5, minor=True)

    axes.grid(which="minor", color="white", linewidth=3.0)

    axes.tick_params(which="both", length=0, pad=10)

    for spine in axes.spines.values():
        spine.set_visible(False)


    for i in range(n_rows):
        for j in range(n_columns):

            label = cell_labels.iat[i, j]

            if not label:
                continue

            # The text has to flip to white once the cell gets dark.
            shade = values[i, j] / vmax

            axes.text(
                j, i,
                label,
                ha="center",
                va="center",
                fontsize=cell_size,
                fontweight="bold",
                color="white" if shade > 0.55 else "#1A1A1A"
            )


    figure.suptitle(
        "Principal components vs. metadata",
        fontsize=title_size,
        fontweight="bold"
    )

    figure.supxlabel(
        "adjusted R² of a one-way ANOVA of the PC scores across each "
        "variable    *  FDR < 0.05     **  FDR < 0.01     ***  FDR < 0.001\n"
        f"categories with fewer than {min_group_size} samples, and samples "
        "with the variable unrecorded, are excluded from that variable's tests",
        fontsize=caption_size,
        color="#555555"
    )

    bar = figure.colorbar(image, ax=axes, shrink=0.9, pad=0.015)

    bar.set_label("adjusted R²", fontsize=legend_size)

    bar.ax.tick_params(labelsize=legend_size * 0.9, length=0)

    bar.outline.set_visible(False)

    return figure


# ---------------------------------------------------------------------


def main():

    args = parse_args()

    variables = args.variables if args.variables else list(DEFAULT_VARIABLES)

    args.outdir.mkdir(parents=True, exist_ok=True)


    eigenvec = read_eigenvec(args.eigenvec)

    available = sum(1 for c in eigenvec.columns if c.startswith("PC"))

    if args.n_pcs > available:
        raise ValueError(
            f"{args.eigenvec.name} holds {available} PCs, "
            f"--n-pcs {args.n_pcs} asks for more"
        )

    pcs = [f"PC{i}" for i in range(1, args.n_pcs + 1)]

    metadata = pd.read_csv(args.metadata)

    data = join_metadata(eigenvec, metadata, variables)


    results = run_tests(data, pcs, variables, args.min_group_size)

    csv_path = args.outdir / f"{args.prefix}.csv"

    # %g rather than %f: the p-values run down to 1e-90 and would all
    # print as 0.000000.
    results.to_csv(csv_path, index=False, float_format="%.10g")


    figure = draw(results, pcs, variables, args.min_group_size, args.font_scale)

    png_path = args.outdir / f"{args.prefix}.png"
    pdf_path = args.outdir / f"{args.prefix}.pdf"

    figure.savefig(png_path, dpi=200, facecolor="white")
    figure.savefig(pdf_path, facecolor="white")

    plt.close(figure)


    print("\nWrote:")
    print(" ", csv_path)
    print(" ", png_path)
    print(" ", pdf_path)


    tested = results.dropna(subset=["association_strength"])

    if tested.empty:
        return

    print("\nStrongest associations:\n")

    strongest = tested.sort_values(
        "association_strength",
        ascending=False
    ).head(10)

    print(
        strongest[
            ["pc", "metadata", "association_strength", "fdr",
             "n_samples", "n_groups"]
        ].to_string(index=False, float_format=lambda v: f"{v:.4g}")
    )


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
