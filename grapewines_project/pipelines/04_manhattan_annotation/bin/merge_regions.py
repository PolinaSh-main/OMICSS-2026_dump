#!/usr/bin/env python3
"""
Merge neighbouring outlier windows into candidate regions.

A single selective sweep usually shows up as a run of adjacent high-FST
windows, not as one isolated window. Treating each window separately
would report the same locus several times, so touching windows (and
windows separated by no more than --max-gap) are collapsed into one
region.

Input:
    <comparison>_outlier_windows.tsv   from window_fst.R

Output:
    <comparison>_candidate_regions.tsv

    region_id  chrom  start  end  length_bp  n_windows  n_snps
    mean_window_fst  max_window_fst
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def parse_args():

    parser = argparse.ArgumentParser(
        description="Merge adjacent outlier windows into regions"
    )

    parser.add_argument(
        "--windows",
        required=True,
        type=Path,
        help="<comparison>_outlier_windows.tsv"
    )

    parser.add_argument(
        "--max-gap",
        type=int,
        default=50_000,
        help="Merge two outlier windows separated by at most this many "
             "bp (default 50000, i.e. one empty window)"
    )

    parser.add_argument(
        "--comparison",
        default=None,
        help="Name used in the output file (default: from --windows)"
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path
    )

    return parser.parse_args()


def chromosome_sort_key(value: str):
    """
    1, 2, ... 10 rather than 1, 10, 2. Non-numeric names sort last,
    alphabetically.
    """

    digits = "".join(c for c in str(value) if c.isdigit())

    return (0, int(digits), "") if digits else (1, 0, str(value))


def merge(windows: pd.DataFrame, max_gap: int) -> pd.DataFrame:

    regions = []

    for chrom, block in windows.groupby("chrom", sort=False):

        block = block.sort_values("window_start")

        start = None
        end = None
        members = []

        def flush():

            if start is None:
                return

            regions.append(
                {
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "length_bp": int(end - start),
                    "n_windows": len(members),
                    "n_snps": int(sum(m["n_snps"] for m in members)),
                    "mean_window_fst": (
                        sum(m["mean_fst"] for m in members) / len(members)
                    ),
                    "max_window_fst": max(m["mean_fst"] for m in members),
                    "max_snp_fst": max(
                        m.get("max_snp_fst", float("nan"))
                        for m in members
                    ),
                }
            )

        for _, window in block.iterrows():

            if start is not None and window["window_start"] - end <= max_gap:

                end = max(end, window["window_end"])
                members.append(window)
                continue

            flush()

            start = window["window_start"]
            end = window["window_end"]
            members = [window]

        flush()

    merged = pd.DataFrame(regions)

    if merged.empty:
        return merged

    merged = (
        merged
        .assign(_chrom_key=merged["chrom"].map(chromosome_sort_key))
        .sort_values(["_chrom_key", "start"])
        .drop(columns="_chrom_key")
        .reset_index(drop=True)
    )


    #
    # A stable, readable identifier. Ranking happens later, so the id
    # must not encode a rank that would change between runs.
    #

    merged.insert(
        0,
        "region_id",
        [
            f"{row.chrom}:{row.start}-{row.end}"
            for row in merged.itertuples()
        ]
    )

    return merged


def main():

    args = parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    comparison = args.comparison or args.windows.name.replace(
        "_outlier_windows.tsv", ""
    )

    windows = pd.read_csv(args.windows, sep="\t")

    required = {"chrom", "window_start", "window_end", "n_snps", "mean_fst"}

    missing = required - set(windows.columns)

    if missing:
        raise ValueError(
            f"{args.windows} is missing columns: {sorted(missing)}"
        )

    windows["chrom"] = windows["chrom"].astype(str)


    regions = merge(windows, args.max_gap)

    output = args.outdir / f"{comparison}_candidate_regions.tsv"

    regions.to_csv(
        output,
        sep="\t",
        index=False,
        float_format="%.6f"
    )

    print(
        f"{len(windows)} outlier windows -> {len(regions)} regions"
    )

    if not regions.empty:
        print()
        print(
            regions
            .sort_values(
                ["max_window_fst", "mean_window_fst"],
                ascending=False
            )
            .head(10)
            [["region_id", "n_windows", "n_snps",
              "mean_window_fst", "max_window_fst"]]
            .to_string(index=False)
        )

    print()
    print(output)


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
