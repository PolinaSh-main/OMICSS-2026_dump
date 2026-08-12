#!/usr/bin/env python3
"""
Create ADMIXTURE barplot for a single Q matrix with metadata tracks.

Fixes vs original
-----------------
- bar align='edge' → fill_between(step='post') for correct alignment
  with metadata imshow AND ~10× faster rendering
- rasterized=True for lightweight PDF output
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


# ------------------------------------------------------------------ #
#  CLI
# ------------------------------------------------------------------ #

def parse_args():

    parser = argparse.ArgumentParser(
        description="Create ADMIXTURE barplot with metadata tracks"
    )

    parser.add_argument("--q",        required=True, type=Path)
    parser.add_argument("--fam",      required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--order",    required=True, type=Path)
    parser.add_argument("--prefix",   required=True)

    return parser.parse_args()


# ------------------------------------------------------------------ #
#  I/O helpers
# ------------------------------------------------------------------ #

def read_fam(path):

    fam = pd.read_csv(path, sep=r"\s+", header=None)

    fam.columns = ["FID", "IID", "PID", "MID", "SEX", "PHENO"]

    return fam


def read_q(path):

    return pd.read_csv(path, sep=r"\s+", header=None)


def find_metadata_id_column(metadata):

    for col in ("Column 1", "IID", "ID", "id", "sample", "Sample"):
        if col in metadata.columns:
            return col

    raise ValueError(
        f"Cannot find sample ID column. "
        f"Available columns: {list(metadata.columns)}"
    )


# ------------------------------------------------------------------ #
#  Data preparation
# ------------------------------------------------------------------ #

def prepare_data(q, fam, metadata, order):

    if len(q) != len(fam):
        raise ValueError(
            f"Q matrix has {len(q)} samples, FAM has {len(fam)}"
        )

    fam["sample"] = fam["IID"].astype(str)
    order["sample"] = order["IID"].astype(str)

    ordered_samples = order["sample"].tolist()

    if not set(ordered_samples).issubset(set(fam["sample"])):
        raise ValueError("Order file contains IDs absent from FAM")

    # Reindex Q by FAM IIDs, then reorder
    q.index = fam["sample"].values
    q = q.loc[ordered_samples]

    # Reindex metadata
    meta_id = find_metadata_id_column(metadata)
    metadata[meta_id] = metadata[meta_id].astype(str)
    metadata = metadata.set_index(meta_id).reindex(ordered_samples)

    return q, metadata


# ------------------------------------------------------------------ #
#  Palette
# ------------------------------------------------------------------ #

def create_palette(n):

    cmap = plt.cm.get_cmap("tab20", max(n, 20))
    return [cmap(i) for i in range(n)]


# ------------------------------------------------------------------ #
#  Plotting
# ------------------------------------------------------------------ #

def plot_admixture(q, metadata, prefix):

    n_samples = len(q)
    k = q.shape[1]

    # ---- figure geometry ----------------------------------------- #

    width = max(12, n_samples / 15)
    n_meta = len(metadata.columns)
    q_height = 5
    meta_height = 0.35

    heights = [q_height] + [meta_height] * n_meta
    total = sum(heights) + 1.0

    fig = plt.figure(figsize=(width, total))

    gs = fig.add_gridspec(
        1 + n_meta, 1,
        height_ratios=heights,
        hspace=0.05
    )

    # ---- ADMIXTURE barplot --------------------------------------- #

    ax = fig.add_subplot(gs[0])

    colors = create_palette(k)

    x = np.arange(n_samples + 1)         # +1 so last bar renders
    bottom = np.zeros(n_samples + 1)

    q_vals = q.values

    for j in range(k):
        col = q_vals[:, j]
        col_ext = np.append(col, 0.0)     # dummy rightmost value
        top = bottom + col_ext

        ax.fill_between(
            x, bottom, top,
            step="post",
            facecolor=colors[j],
            edgecolor="none",
            linewidth=0,
            rasterized=True
        )

        bottom = top

    ax.set_xlim(0, n_samples)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Ancestry proportion")
    ax.set_title(f"ADMIXTURE K={k}")
    ax.set_xticks([])

    handles = [
        plt.Rectangle((0, 0), 1, 1, fc=colors[i])
        for i in range(k)
    ]
    labels = [f"K{i+1}" for i in range(k)]

    ax.legend(
        handles, labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.15),
        ncol=k,
        frameon=False,
        fontsize=8
    )

    # ---- metadata tracks ----------------------------------------- #

    for row, column in enumerate(metadata.columns, start=1):

        meta_ax = fig.add_subplot(gs[row], sharex=ax)

        values = metadata[column].fillna("NA").astype(str)
        categories = values.unique().tolist()
        mapping = {v: i for i, v in enumerate(categories)}
        numeric = np.array([mapping[v] for v in values])

        cmap = ListedColormap(create_palette(len(categories)))

        meta_ax.imshow(
            [numeric],
            aspect="auto",
            interpolation="nearest",
            cmap=cmap,
            extent=[0, n_samples, 0, 1],
            rasterized=True
        )

        meta_ax.set_ylabel(
            column,
            rotation=0,
            ha="right",
            va="center",
            fontsize=8
        )
        meta_ax.yaxis.set_label_coords(-0.02, 0.5)
        meta_ax.set_yticks([])
        meta_ax.set_xticks([])

    # ---- save ---------------------------------------------------- #

    plt.subplots_adjust(left=0.22, right=0.98, top=0.95, bottom=0.05)

    fig.savefig(f"{prefix}.png", dpi=200, bbox_inches="tight")
    fig.savefig(f"{prefix}.pdf",           bbox_inches="tight")

    plt.close(fig)


# ------------------------------------------------------------------ #
#  Entry
# ------------------------------------------------------------------ #

def main():

    args = parse_args()

    q        = read_q(args.q)
    fam      = read_fam(args.fam)
    metadata = pd.read_csv(args.metadata)
    order    = pd.read_csv(args.order, sep="\t")

    q, metadata = prepare_data(q, fam, metadata, order)

    plot_admixture(q, metadata, args.prefix)


if __name__ == "__main__":
    main()