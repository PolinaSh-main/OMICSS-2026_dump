#!/usr/bin/env python3

import argparse
from pathlib import Path

import scanpy as sc
import anndata as ad
import gzip


def read_10x_custom(input_dir: Path):

    """
    Reads 10x matrix with 2-column features.tsv:
    gene_id <tab> gene_symbol
    """

    matrix = input_dir / "matrix.mtx.gz"
    barcodes = input_dir / "barcodes.tsv.gz"
    features = input_dir / "features.tsv.gz"

    if not matrix.exists():
        raise FileNotFoundError(matrix)

    if not barcodes.exists():
        raise FileNotFoundError(barcodes)

    if not features.exists():
        raise FileNotFoundError(features)

    adata = sc.read_mtx(matrix).T

    # barcodes
    with gzip.open(barcodes, "rt") as f:
        adata.obs_names = [
            x.strip()
            for x in f
        ]

    # features
    genes = []

    with gzip.open(features, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            genes.append(parts)

    gene_ids = [x[0] for x in genes]
    gene_symbols = [x[1] for x in genes]

    adata.var_names = gene_symbols
    adata.var["gene_ids"] = gene_ids
    adata.var["gene_symbols"] = gene_symbols

    return adata



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
        "--condition",
        required=True
    )

    parser.add_argument(
        "--min-genes",
        type=int,
        default=200
    )

    parser.add_argument(
        "--max-genes",
        type=int,
        default=2500
    )

    parser.add_argument(
        "--max-mt",
        type=float,
        default=5.0
    )


    args = parser.parse_args()


    input_dir = Path(args.input)


    print(f"Reading {input_dir}")

    adata = read_10x_custom(input_dir)


    adata.obs["condition"] = args.condition


    # QC metrics

    adata.var["mt"] = (
        adata.var_names
        .str.upper()
        .str.startswith("MT-")
    )


    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        inplace=True
    )


    print(
        "Before filtering:",
        adata.n_obs,
        "cells"
    )


    # filtering

    adata = adata[
        (adata.obs.n_genes_by_counts >= args.min_genes)
        &
        (adata.obs.n_genes_by_counts <= args.max_genes)
        &
        (adata.obs.pct_counts_mt <= args.max_mt)
    ].copy()


    print(
        "After filtering:",
        adata.n_obs,
        "cells"
    )


    output = Path(args.output) / "qc_filtered.h5ad"

    adata.write_h5ad(output)


if __name__ == "__main__":
    main()