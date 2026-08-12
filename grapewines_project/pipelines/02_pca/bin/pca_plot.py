#!/usr/bin/env python3

import argparse
import os
import re

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


sns.set_theme(style="whitegrid")


def clean_name(name):
    return re.sub(r"[^A-Za-z0-9_]+", "_", name)


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--eigenvec", required=True)
    parser.add_argument("--eigenval", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--color-column", required=True)
    parser.add_argument("--outdir", required=True)

    args = parser.parse_args()


    os.makedirs(args.outdir, exist_ok=True)


    # -----------------------
    # Read PCA coordinates
    # -----------------------

    pcs = pd.read_csv(
        args.eigenvec,
        sep=r"\s+",
        header=None
    )

    n_pc = pcs.shape[1] - 2

    pcs.columns = (
        ["FID", "IID"] +
        [f"PC{i}" for i in range(1, n_pc + 1)]
    )


    # -----------------------
    # Eigenvalues
    # -----------------------

    eigenvalues = pd.read_csv(
        args.eigenval,
        header=None
    )[0]

    explained = eigenvalues / eigenvalues.sum()


    # -----------------------
    # Metadata
    # -----------------------

    meta = pd.read_csv(args.metadata)

    meta = meta.rename(
        columns={
            "Column 1": "IID"
        }
    )


    if args.color_column not in meta.columns:
        raise ValueError(
            f"Column '{args.color_column}' not found.\n"
            f"Available columns:\n{list(meta.columns)}"
        )


    pcs = pcs.merge(
        meta[["IID", args.color_column]],
        on="IID",
        how="left"
    )


    # -----------------------
    # Plot function
    # -----------------------

    def plot_pc(x, y):

        plt.figure(figsize=(8,6))

        sns.scatterplot(
            data=pcs,
            x=x,
            y=y,
            hue=args.color_column,
            s=60
        )


        plt.xlabel(
            f"{x} ({explained[int(x[2:])-1]*100:.2f}%)"
        )

        plt.ylabel(
            f"{y} ({explained[int(y[2:])-1]*100:.2f}%)"
        )


        plt.title(
            f"{x} vs {y}\ncolored by {args.color_column}"
        )


        plt.legend(
            bbox_to_anchor=(1.05,1),
            loc="upper left"
        )


        plt.tight_layout()


        plt.savefig(
            os.path.join(
                args.outdir,
                f"{x}_{y}.png"
            ),
            dpi=300,
            bbox_inches="tight"
        )

        plt.close()


    plot_pc("PC1","PC2")
    plot_pc("PC1","PC3")
    plot_pc("PC2","PC3")


if __name__ == "__main__":
    main()