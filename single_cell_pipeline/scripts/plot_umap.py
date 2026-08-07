#!/usr/bin/env python3

from pathlib import Path

import scanpy as sc


adata = sc.read_h5ad(
    "results/umap/umap.h5ad"
)


out = Path("results/figures")
out.mkdir(
    parents=True,
    exist_ok=True
)


sc.settings.figdir = out


# condition
sc.pl.umap(
    adata,
    color="condition",
    save="_condition.png",
    show=False
)


# clusters
sc.pl.umap(
    adata,
    color="leiden",
    save="_clusters.png",
    show=False
)


# QC
for col in [
    "n_genes",
    "total_counts",
    "pct_counts_mt"
]:
    if col in adata.obs:
        sc.pl.umap(
            adata,
            color=col,
            save=f"_{col}.png",
            show=False
        )