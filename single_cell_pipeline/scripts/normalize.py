#!/usr/bin/env python3

import argparse
from pathlib import Path

import scanpy as sc
import anndata as ad


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        required=True
    )

    parser.add_argument(
        "--output",
        required=True
    )

    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=2000
    )

    parser.add_argument(
        "--n-pcs",
        type=int,
        default=30
    )

    args = parser.parse_args()


    adata = ad.read_h5ad(args.input)


    print("Normalizing")


    sc.pp.normalize_total(
        adata,
        target_sum=1e4
    )

    sc.pp.log1p(adata)


    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=args.n_top_genes,
        flavor="seurat_v3"
    )


    adata = adata[
        :,
        adata.var.highly_variable
    ].copy()


    sc.pp.scale(
        adata,
        max_value=10
    )


    sc.tl.pca(
        adata,
        n_comps=args.n_pcs
    )


    adata.write_h5ad(
        Path(args.output) / "normalized_pca.h5ad"
    )


if __name__ == "__main__":
    main()