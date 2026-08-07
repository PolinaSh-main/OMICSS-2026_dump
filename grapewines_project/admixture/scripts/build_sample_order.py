#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fam",
        required=True,
        type=Path
    )

    parser.add_argument(
        "--q",
        required=True,
        type=Path
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path
    )

    return parser.parse_args()


def main():

    args = parse_args()

    #
    # Read PLINK fam
    #
    fam = pd.read_csv(
        args.fam,
        sep=r"\s+",
        header=None
    )

    fam.columns = [
        "FID",
        "IID",
        "PID",
        "MID",
        "SEX",
        "PHENO"
    ]

    #
    # Read K=2 Q matrix
    #
    q = pd.read_csv(
        args.q,
        sep=r"\s+",
        header=None
    )

    if len(fam) != len(q):

        raise ValueError(
            f"Number of samples differs: "
            f"fam={len(fam)}, Q={len(q)}"
        )


    #
    # First ancestry component
    #
    fam["K1"] = q.iloc[:, 0]


    #
    # Sort by K1
    #
    fam = fam.sort_values(
        "K1"
    )


    #
    # Save order
    #
    order = fam[
        [
            "FID",
            "IID"
        ]
    ]


    order.to_csv(
        args.output,
        sep="\t",
        index=False
    )


    print(
        f"Saved sample order for {len(order)} samples"
    )


if __name__ == "__main__":

    main()