#!/usr/bin/env python3
"""
Assign samples to ADMIXTURE groups.

A sample joins group K<i> when its largest ancestry proportion is in
column i AND that proportion is at least --min-q. Everything else is
admixed and goes to a separate file.

Input:
    .fam   sample IDs, in the same row order as the Q file
    .Q     ADMIXTURE ancestry proportions, one row per sample

Output:
    K<i>_samples.txt          one sample ID per line, for vcftools
    admixed_samples.txt       samples below the threshold
    sample_q_values_K<K>.tsv  every sample with its Q values and assignment
    group_sizes_K<K>.tsv      how many samples ended up in each group
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def parse_args():

    parser = argparse.ArgumentParser(
        description="Assign samples to ADMIXTURE groups by max Q"
    )

    parser.add_argument(
        "--fam",
        required=True,
        type=Path,
        help="PLINK .fam, row order must match the Q file"
    )

    parser.add_argument(
        "--q",
        required=True,
        type=Path,
        help="ADMIXTURE .Q file"
    )

    parser.add_argument(
        "--min-q",
        type=float,
        default=0.75,
        help="Minimum ancestry proportion to assign a sample (default 0.75)"
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path
    )

    return parser.parse_args()


def read_fam(path: Path) -> pd.DataFrame:

    fam = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        dtype=str
    )

    if fam.shape[1] < 2:
        raise ValueError(
            f"{path} does not look like a .fam file "
            f"({fam.shape[1]} columns)"
        )

    fam = fam.iloc[:, :2]

    fam.columns = [
        "FID",
        "IID"
    ]

    return fam


def read_q(path: Path) -> pd.DataFrame:

    q = pd.read_csv(
        path,
        sep=r"\s+",
        header=None
    )

    q.columns = [
        f"K{i + 1}"
        for i in range(q.shape[1])
    ]

    return q


def main():

    args = parse_args()

    args.outdir.mkdir(
        parents=True,
        exist_ok=True
    )


    fam = read_fam(args.fam)

    q = read_q(args.q)


    #
    # The Q file carries no sample names at all -- it is matched to the
    # .fam purely by row order. If the two disagree every downstream
    # result would be silently wrong, so refuse to continue.
    #

    if len(fam) != len(q):

        raise ValueError(
            f"Row count mismatch: "
            f"{args.fam.name} has {len(fam)} samples, "
            f"{args.q.name} has {len(q)} rows"
        )


    k = q.shape[1]

    group_cols = list(q.columns)


    df = pd.concat(
        [
            fam.reset_index(drop=True),
            q.reset_index(drop=True)
        ],
        axis=1
    )


    df["max_q"] = df[group_cols].max(axis=1)

    df["best_group"] = df[group_cols].idxmax(axis=1)


    df["assignment"] = df["best_group"].where(
        df["max_q"] >= args.min_q,
        other="admixed"
    )


    #
    # Sample lists, one file per group. vcftools wants bare sample IDs.
    #

    for group in group_cols:

        members = df.loc[
            df["assignment"] == group,
            "IID"
        ]

        out = args.outdir / f"{group}_samples.txt"

        members.to_csv(
            out,
            index=False,
            header=False
        )

        print(f"{group:>4}  {len(members):4d} samples  -> {out.name}")


    admixed = df.loc[
        df["assignment"] == "admixed",
        "IID"
    ]

    admixed.to_csv(
        args.outdir / "admixed_samples.txt",
        index=False,
        header=False
    )

    print(f"{'adm':>4}  {len(admixed):4d} samples  -> admixed_samples.txt")


    #
    # Per-sample table: the audit trail for the assignment above, and
    # the input for naming the groups from metadata.
    #

    df.drop(columns=["best_group"]).to_csv(
        args.outdir / f"sample_q_values_K{k}.tsv",
        sep="\t",
        index=False,
        float_format="%.6f"
    )


    #
    # Group sizes -- deliverable in its own right.
    #

    sizes = (
        df["assignment"]
        .value_counts()
        .reindex(group_cols + ["admixed"])
        .fillna(0)
        .astype(int)
        .rename_axis("group")
        .reset_index(name="n_samples")
    )

    sizes["comment"] = [
        "admixed, excluded from FST" if g == "admixed"
        else "small group, interpret with care" if n < 5
        else "ok"
        for g, n in zip(sizes["group"], sizes["n_samples"])
    ]

    sizes.to_csv(
        args.outdir / f"group_sizes_K{k}.tsv",
        sep="\t",
        index=False
    )

    print()
    print(sizes.to_string(index=False))


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
