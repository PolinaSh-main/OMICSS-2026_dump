#!/usr/bin/env python3
"""
Plot ADMIXTURE cross-validation error curve.
"""

import argparse

import pandas as pd
import matplotlib.pyplot as plt


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        required=True,
        help="TSV file produced by extract_cv.py"
    )

    parser.add_argument(
        "--png",
        required=True,
        help="Output PNG"
    )

    parser.add_argument(
        "--pdf",
        required=True,
        help="Output PDF"
    )

    return parser.parse_args()


def main():

    args = parse_args()

    df = pd.read_csv(args.input, sep="\t")

    best = df.loc[df["CV"].idxmin()]

    fig, ax = plt.subplots(figsize=(6, 4))

    ax.plot(
        df["K"].to_numpy(),
        df["CV"].to_numpy(),
        marker="o",
        linewidth=2,
        zorder=2
    )

    ax.scatter(
        [best["K"]],
        [best["CV"]],
        s=120,
        color="tab:red",
        zorder=3
    )

    ax.annotate(
        f'Best K={int(best["K"])}',
        (best["K"], best["CV"]),
        xytext=(10, 10),
        textcoords="offset points",
        fontweight="bold"
    )

    ax.set_xlabel("Number of Ancestral Populations (K)")
    ax.set_ylabel("Cross-Validation (CV) Error")
    ax.set_title("Cross-Validation Error Across Different K Values")

    ax.set_xticks(df["K"])

    ax.grid(True, alpha=0.4)

    fig.tight_layout()

    fig.savefig(args.png, dpi=300)
    fig.savefig(args.pdf)

    plt.close(fig)


if __name__ == "__main__":
    main()