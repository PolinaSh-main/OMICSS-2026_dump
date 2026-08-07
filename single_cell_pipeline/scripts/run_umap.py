#!/usr/bin/env python3

import argparse
from pathlib import Path

import scanpy as sc


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

    args = parser.parse_args()


    output = Path(args.output)
    output.mkdir(
        parents=True,
        exist_ok=True
    )


    adata = sc.read_h5ad(args.input)


    print("Building neighbors")

    sc.pp.neighbors(
        adata,
        use_rep="X_harmony",
        n_neighbors=15,
        n_pcs=30
    )


    print("Calculating UMAP")

    sc.tl.umap(
        adata
    )


    print("Clustering")

    sc.tl.leiden(
        adata,
        resolution=0.5
    )


    adata.write_h5ad(
        output / "umap.h5ad"
    )


if __name__ == "__main__":
    main()