#!/usr/bin/env python3
"""
PC1/PC2 (and PC1/PC3) coloured by the ADMIXTURE K = 7 assignment.

The by_metadata plots colour the same coordinates by passport data --
country, subspecies, berry colour. This one colours them by the groups
we built ourselves at Q >= 0.75, so PCA and ADMIXTURE can be read
against each other on one picture instead of two.

Admixed samples are drawn first, in grey, so the 160 assigned ones stay
visible on top of the 252 that are not.

Colours match the tree figures in 05_tree.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd


K_COLOURS = {
    "K1": "#1F77B4",
    "K2": "#FF7F0E",
    "K3": "#2CA02C",
    "K4": "#BCBD22",
    "K5": "#D62728",
    "K6": "#9467BD",
    "K7": "#17BECF",
}

ADMIXED_COLOUR = "#BDBDBD"


def parse_args():

    here = Path(__file__).resolve().parents[3]

    parser = argparse.ArgumentParser(
        description="PCA scatter coloured by ADMIXTURE assignment"
    )

    parser.add_argument(
        "--coordinates",
        type=Path,
        default=here / "results/pca/global/pca_coordinates.tsv"
    )

    parser.add_argument(
        "--assignments",
        type=Path,
        default=here / "results/fst/K7/sample_lists/sample_q_values_K7.tsv"
    )

    parser.add_argument(
        "--variance",
        type=Path,
        default=here / "results/pca/global/explained_variance.tsv"
    )

    parser.add_argument(
        "--outdir",
        type=Path,
        default=here / "results/pca/by_admixture"
    )

    return parser.parse_args()


def load(args) -> tuple[pd.DataFrame, dict]:

    coords = pd.read_csv(args.coordinates, sep="\t")

    assign = pd.read_csv(
        args.assignments,
        sep="\t",
        usecols=["IID", "max_q", "assignment"]
    )

    merged = coords.merge(assign, on="IID", how="left", validate="1:1")

    missing = int(merged["assignment"].isna().sum())

    if missing:
        raise ValueError(
            f"{missing} samples in the PCA table have no ADMIXTURE "
            f"assignment; the two runs are not on the same sample set"
        )

    variance = pd.read_csv(args.variance, sep="\t")

    percent = {
        row.PC: 100 * row.Explained_variance
        for row in variance.itertuples()
    }

    return merged, percent


def draw(ax, data: pd.DataFrame, x: str, y: str, percent: dict):

    admixed = data[data["assignment"] == "admixed"]

    ax.scatter(
        admixed[x],
        admixed[y],
        s=26,
        c=ADMIXED_COLOUR,
        edgecolors="white",
        linewidths=0.4,
        label=f"Admixed (n = {len(admixed)})",
        zorder=1
    )

    for group in sorted(K_COLOURS):

        subset = data[data["assignment"] == group]

        if subset.empty:
            continue

        ax.scatter(
            subset[x],
            subset[y],
            s=52,
            c=K_COLOURS[group],
            edgecolors="black",
            linewidths=0.5,
            label=f"{group} (n = {len(subset)})",
            zorder=2
        )

    ax.set_xlabel(f"{x} ({percent.get(x, float('nan')):.2f}%)")
    ax.set_ylabel(f"{y} ({percent.get(y, float('nan')):.2f}%)")
    ax.grid(True, linewidth=0.4, alpha=0.5)
    ax.set_axisbelow(True)


def main():

    args = parse_args()

    data, percent = load(args)

    args.outdir.mkdir(parents=True, exist_ok=True)

    for x, y in [("PC1", "PC2"), ("PC1", "PC3"), ("PC2", "PC3")]:

        figure, ax = plt.subplots(figsize=(9.5, 7.5))

        draw(ax, data, x, y, percent)

        ax.set_title(
            f"{x} vs {y}, coloured by ADMIXTURE K = 7 "
            f"(assigned at Q ≥ 0.75)"
        )

        ax.legend(
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=True,
            fontsize=9
        )

        figure.tight_layout()

        for suffix in ("png", "pdf"):
            figure.savefig(
                args.outdir / f"{x}_{y}.{suffix}",
                dpi=150,
                bbox_inches="tight"
            )

        plt.close(figure)

        print(f"written {args.outdir / f'{x}_{y}.png'}")


if __name__ == "__main__":
    main()
