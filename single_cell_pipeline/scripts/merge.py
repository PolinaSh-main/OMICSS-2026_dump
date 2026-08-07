#!/usr/bin/env python3

import argparse
from pathlib import Path

import anndata as ad


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True
    )

    parser.add_argument(
        "--output",
        required=True
    )

    args = parser.parse_args()


    objects = []

    for file in args.inputs:

        path = Path(file)

        print(f"Reading {path}")

        adata = ad.read_h5ad(path)
        
        adata.var_names_make_unique()

        objects.append(adata)


    print("Merging")

    merged = ad.concat(
        objects,
        join="outer",
        label="batch",
        keys=[
            x.obs["condition"].iloc[0]
            if "condition" in x.obs
            else f"sample_{i}"
            for i, x in enumerate(objects)
        ],
        index_unique="-"
    )


    output = Path(args.output)

    merged.write_h5ad(output)

    print(
        "Merged:",
        merged.shape
    )


if __name__ == "__main__":
    main()