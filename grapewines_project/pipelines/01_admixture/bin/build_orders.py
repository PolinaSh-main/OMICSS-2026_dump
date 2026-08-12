#!/usr/bin/env python3

import argparse
from pathlib import Path
import re

import pandas as pd
import numpy as np



def parse_args():

    p = argparse.ArgumentParser()

    p.add_argument(
        "--qdir",
        type=Path,
        required=True
    )

    p.add_argument(
        "--fam",
        type=Path,
        required=True
    )

    p.add_argument(
        "--metadata",
        type=Path,
        required=True
    )

    p.add_argument(
        "--outdir",
        type=Path,
        required=True
    )

    return p.parse_args()



def get_k(path):

    m = re.search(
        r"\.final\.(\d+)\.Q",
        path.name
    )

    return int(m.group(1))



def read_q(path):

    return pd.read_csv(
        path,
        sep=r"\s+",
        header=None
    )



def read_fam(path):

    fam = pd.read_csv(
        path,
        sep=r"\s+",
        header=None
    )

    fam.columns = [
        "FID",
        "IID",
        "father",
        "mother",
        "sex",
        "phenotype"
    ]

    return fam[["FID", "IID"]]

def find_country_column(metadata):

    for c in metadata.columns:

        if "country" in c.lower():

            return c

    raise ValueError(
        "Country column not found"
    )



def main():

    args = parse_args()

    args.outdir.mkdir(
        exist_ok=True,
        parents=True
    )


    fam = read_fam(
        args.fam
    )


    metadata = pd.read_csv(
        args.metadata
    )


    country_col = find_country_column(
        metadata
    )


    #
    # Merge metadata with sample order
    #

    meta = fam.merge(
        metadata,
        left_on="IID",
        right_on="Column 1",
        how="left"
    )


    if meta[country_col].isna().all():

        raise ValueError(
            "Metadata IDs do not match FAM"
        )


    q_files = sorted(
        args.qdir.glob(
            "cauc_filtered.final.*.Q"
        ),
        key=get_k
    )


    for q_file in q_files:


        k = get_k(
            q_file
        )


        q = read_q(
            q_file
        )


        if len(q) != len(meta):

            raise ValueError(
                f"K={k}: "
                f"Q rows {len(q)} != samples {len(meta)}"
            )


        q.columns = [
            f"Q{i+1}"
            for i in range(k)
        ]


        df = pd.concat(
            [
                meta,
                q
            ],
            axis=1
        )


        #
        # Major ancestry component
        #

        q_cols = [
            f"Q{i+1}"
            for i in range(k)
        ]


        df["major_Q"] = (
            df[q_cols]
            .idxmax(axis=1)
        )


        df["major_value"] = (
            df[q_cols]
            .max(axis=1)
        )


        #
        # Sort:
        #
        # 1. major Q cluster
        # 2. country
        # 3. decreasing major component
        #

        df = df.sort_values(
            [
                "major_Q",
                country_col,
                "major_value"
            ],
            ascending=[
                True,
                True,
                False
            ]
        )


        df[
            [
                "FID",
                "IID"
            ]
        ].to_csv(
            args.outdir /
            f"order_Q{k}.tsv",
            sep="\t",
            index=False
        )



if __name__ == "__main__":
    main()