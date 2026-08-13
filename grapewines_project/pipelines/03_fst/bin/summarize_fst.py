#!/usr/bin/env python3
"""
Collapse per-SNP FST files into one table, one row per comparison.

Each *.weir.fst holds several hundred thousand rows -- one per SNP --
which is far too much to read directly. This reduces every comparison
to a handful of numbers.

Two averages are reported, because they answer different questions:

    mean_fst      unweighted average over SNPs. This is the number the
                  course instructions ask for.
    weighted_fst  the ratio-of-averages estimator vcftools prints in
                  its log. It is less sensitive to the many
                  low-information SNPs and is usually the one quoted in
                  papers.

Sites where vcftools could not compute FST (monomorphic within both
groups) come out as -nan and are excluded from the mean.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


COMPARISON_RE = re.compile(r"^(K\d+)_vs_(K\d+)$")

MEAN_RE = re.compile(
    r"Weir and Cockerham mean Fst estimate:\s*([0-9.eE+-]+)"
)

WEIGHTED_RE = re.compile(
    r"Weir and Cockerham weighted Fst estimate:\s*([0-9.eE+-]+)"
)


def parse_args():

    parser = argparse.ArgumentParser(
        description="Summarise pairwise FST results"
    )

    parser.add_argument(
        "--fstdir",
        required=True,
        type=Path,
        help="Directory with *.weir.fst (and optional vcftools *.log)"
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path
    )

    return parser.parse_args()


def summarise_one(path: Path) -> dict:

    name = path.name.replace(".weir.fst", "")

    match = COMPARISON_RE.match(name)

    if match is None:
        raise ValueError(
            f"Cannot read group names from {path.name}; "
            f"expected K<i>_vs_K<j>.weir.fst"
        )

    group_a, group_b = match.groups()


    fst = pd.read_csv(
        path,
        sep="\t",
        usecols=["WEIR_AND_COCKERHAM_FST"],
        na_values=["-nan", "nan", "NaN", "-inf", "inf"]
    )["WEIR_AND_COCKERHAM_FST"]

    usable = fst.dropna()


    row = {
        "comparison": name,
        "group_a": group_a,
        "group_b": group_b,
        "n_sites": int(len(fst)),
        "n_sites_used": int(len(usable)),
        "mean_fst": float(usable.mean()) if len(usable) else np.nan,
        #
        # Weir and Cockerham FST goes negative wherever the within-group
        # variance happens to exceed the between-group variance, which
        # at a single site is noise around zero rather than a real
        # quantity. Two conventions are in use: average the estimates as
        # they come (mean_fst -- this is also what vcftools prints in
        # its own log), or set the negatives to zero first
        # (mean_fst_nonneg). The second is always the larger of the two,
        # by about 3-4% on this data, because it truncates one tail of
        # the noise and not the other.
        #
        # Both are reported so that a figure can say which it used, and
        # so that two figures built on different conventions cannot
        # quietly disagree in a slide deck.
        #
        "mean_fst_nonneg": (
            float(usable.clip(lower=0).mean()) if len(usable) else np.nan
        ),
        "median_fst": float(usable.median()) if len(usable) else np.nan,
        "max_fst": float(usable.max()) if len(usable) else np.nan,
        "weighted_fst": np.nan,
    }


    #
    # vcftools already computed both estimates; take the weighted one
    # from its log rather than reimplementing the estimator.
    #

    log = path.with_name(f"{name}.log")

    if log.exists():

        text = log.read_text(errors="ignore")

        weighted = WEIGHTED_RE.search(text)

        if weighted:
            row["weighted_fst"] = float(weighted.group(1))

        mean = MEAN_RE.search(text)

        if mean:
            #
            # Sanity check: our own mean should agree with vcftools.
            #
            reported = float(mean.group(1))

            if (
                not np.isnan(row["mean_fst"])
                and abs(reported - row["mean_fst"]) > 1e-3
            ):
                print(
                    f"WARNING: {name}: mean FST {row['mean_fst']:.5f} "
                    f"differs from the vcftools log ({reported:.5f})",
                    file=sys.stderr
                )

    return row


def main():

    args = parse_args()

    files = sorted(args.fstdir.glob("*.weir.fst"))

    if not files:
        raise FileNotFoundError(
            f"No *.weir.fst files in {args.fstdir}"
        )


    rows = [
        summarise_one(path)
        for path in files
    ]


    summary = pd.DataFrame(rows)


    def sort_key(column):
        return summary[column].str.slice(1).astype(int)

    summary = (
        summary
        .assign(
            _a=sort_key("group_a"),
            _b=sort_key("group_b")
        )
        .sort_values(["_a", "_b"])
        .drop(columns=["_a", "_b"])
        .reset_index(drop=True)
    )


    args.output.parent.mkdir(parents=True, exist_ok=True)

    summary.to_csv(
        args.output,
        sep="\t",
        index=False,
        float_format="%.6f"
    )


    print(
        summary[
            ["comparison", "n_sites_used", "mean_fst", "weighted_fst"]
        ].to_string(index=False)
    )

    print()
    print(f"Written {len(summary)} comparisons to {args.output}")


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
