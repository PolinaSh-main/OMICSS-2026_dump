#!/usr/bin/env python3
"""
Draw equal-sized random subsets of every ADMIXTURE group.

The pairwise FST ranking from 03_fst tracks group size (Spearman -0.79)
rather than anything biological. FST is between-group variance divided
by total variance, and the groups are not comparable objects: at
Q >= 0.75 a component with 12 members contributes its tight core, while
one with 49 contributes a broad slice. A tight core has little variance
inside it, so the denominator shrinks and FST rises.

Cutting every group to the size of the smallest makes the sample sizes
comparable. It does NOT make the groups equally diverse -- twelve
samples drawn from a broad group of 49 stay broader than the twelve that
are all of a narrow group. So this is a diagnostic, not a correction:

    the correlation with size disappears  -> it was a sampling effect
                                             and the new ranking means
                                             something
    the correlation survives              -> the narrow groups really
                                             are narrow, and the cause
                                             is the Q threshold rather
                                             than the sample size

One draw is itself arbitrary, so several independent draws are made and
the spread between them is reported alongside the mean.
"""

from __future__ import annotations

import argparse
import zlib
from pathlib import Path

import pandas as pd


def parse_args():

    parser = argparse.ArgumentParser(
        description="Equal-size random subsets of each ancestry group"
    )

    parser.add_argument("--assignments", required=True, type=Path)

    parser.add_argument(
        "--n",
        type=int,
        default=0,
        help="Samples per group; 0 means use the smallest group"
    )

    parser.add_argument("--replicates", type=int, default=10)

    parser.add_argument("--seed", type=int, default=20260813)

    parser.add_argument("--outdir", type=Path, default=Path("."))

    return parser.parse_args()


def main():

    args = parse_args()

    table = pd.read_csv(args.assignments, sep="\t")

    assigned = table[table["assignment"] != "admixed"]

    sizes = assigned["assignment"].value_counts().sort_index()

    if sizes.empty:
        raise ValueError("No assigned samples; nothing to subsample")

    n = args.n if args.n > 0 else int(sizes.min())

    too_small = sizes[sizes < n]

    if not too_small.empty:
        raise ValueError(
            f"n = {n} exceeds the size of {dict(too_small)}"
        )

    args.outdir.mkdir(parents=True, exist_ok=True)


    #
    # The full assigned list, used once to cut the VCF down from 412
    # samples to 160 so that the 21 x replicates vcftools passes each
    # read a much smaller file.
    #

    (args.outdir / "assigned_samples.txt").write_text(
        "\n".join(assigned["IID"].astype(str)) + "\n",
        encoding="utf-8"
    )

    print(f"{len(assigned)} assigned samples in {len(sizes)} groups")
    print(f"drawing {n} per group, {args.replicates} replicates\n")
    print(sizes.to_string())
    print()


    for replicate in range(1, args.replicates + 1):

        for group, members in assigned.groupby("assignment"):

            #
            # Seed per (replicate, group) so a rerun reproduces exactly
            # the same draw, and so adding a replicate does not shift
            # the ones already computed.
            #
            # crc32 rather than hash(): Python randomises string hashing
            # per process unless PYTHONHASHSEED is set, which would make
            # every rerun draw a different subset.
            #
            drawn = members.sample(
                n=n,
                random_state=(
                    args.seed
                    + 1000 * replicate
                    + zlib.crc32(group.encode()) % 997
                )
            )

            path = args.outdir / f"rep{replicate:02d}_{group}_samples.txt"

            path.write_text(
                "\n".join(drawn["IID"].astype(str)) + "\n",
                encoding="utf-8"
            )

    print(
        f"written {args.replicates * len(sizes)} sample lists "
        f"to {args.outdir}"
    )


if __name__ == "__main__":
    main()
