#!/usr/bin/env python3
"""
Two sanity checks on the windowed-FST outliers.

Neither is a proper neutrality test -- that needs a demographic model we
do not have. Both attack the same weakness from the side, using only
tables that already exist.


CHECK 1 -- does an outlier window recur across group pairs?

Each of the 21 comparisons keeps its own top windows. A window driven by
selection in, say, the cultivated lineage should surface in several
comparisons that involve cultivated groups. A window that is high
because one 12-sample group happened to drift there should surface once
and never again.

The null here is cheap and honest, because it does not need to model
evolution at all: if each pair marked its windows at random, how often
would the same window be marked by k pairs? Pair p marks W_p windows out
of N; scatter those at random, repeat, count. Anything far above that
line is real recurrence.

Recurrence alone is not enough, though. If a window is an outlier in
K1_vs_K2, K1_vs_K3, K1_vs_K4 ... that is not six independent
observations -- it is one statement about K1. So recurrent windows are
split into

    group-driven   every pair that flags it shares one group
    shared         flagged by pairs with no group in common, i.e. at
                   least two disjoint comparisons

Only the second kind is evidence of something bigger than one group.


CHECK 2 -- is the window mean driven by how few SNPs are in it?

A 50 kb window holding 6 SNPs has a mean FST with a huge standard error;
one holding 60 has a stable one. If the top of the distribution is
populated by SNP-poor windows, the ranking is measuring sampling noise
rather than differentiation.

Without the full genome-wide window table we cannot compare outliers to
background, so the test used here is internal: among the outliers, does
window FST go up as SNP count goes down? A strong negative correlation
is the signature of the artefact.
"""

from __future__ import annotations

import argparse
import itertools
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


WINDOW_BP = 50_000

PAIR_RE = re.compile(r"^(K\d+)_vs_(K\d+)$")


def parse_args():

    here = Path(__file__).resolve().parents[3]

    parser = argparse.ArgumentParser(
        description="Recurrence and SNP-density checks on FST outliers"
    )

    parser.add_argument(
        "--manhattan-dir",
        type=Path,
        default=here / "results/manhattan/K7"
    )

    parser.add_argument(
        "--outdir",
        type=Path,
        default=here / "results/manhattan/K7/_checks"
    )

    parser.add_argument(
        "--total-windows",
        type=int,
        default=8387,
        help="Genome-wide 50 kb windows that carried enough SNPs to be "
             "scored; the denominator for the random-scatter null"
    )

    parser.add_argument(
        "--replicates",
        type=int,
        default=2000
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=20260813
    )

    return parser.parse_args()


def load_regions(manhattan_dir: Path) -> pd.DataFrame:
    """One row per (pair, 50 kb window)."""

    rows = []

    for path in sorted(manhattan_dir.glob("*/*_candidate_regions.tsv")):

        pair = path.name.replace("_candidate_regions.tsv", "")

        if PAIR_RE.match(pair) is None:
            continue

        table = pd.read_csv(path, sep="\t")

        for region in table.itertuples():

            #
            # Regions are runs of adjacent outlier windows, merged. The
            # merged boundaries differ between pairs, so matching on
            # region_id would miss real overlaps. Split back into the
            # fixed 50 kb grid and match on that.
            #
            for start in range(
                int(region.start), int(region.end), WINDOW_BP
            ):
                rows.append({
                    "pair": pair,
                    "group_a": PAIR_RE.match(pair).group(1),
                    "group_b": PAIR_RE.match(pair).group(2),
                    "window": f"{region.chrom}:{start}",
                    "chrom": int(region.chrom),
                    "start": start,
                    "n_snps": int(region.n_snps),
                    "n_windows": int(region.n_windows),
                    "mean_window_fst": float(region.mean_window_fst),
                    "max_window_fst": float(region.max_window_fst),
                })

    if not rows:
        raise FileNotFoundError(
            f"No *_candidate_regions.tsv under {manhattan_dir}"
        )

    return pd.DataFrame(rows)


def classify(pairs: list[str]) -> str:
    """
    'group-driven' if some single group appears in every pair that
    flagged this window; 'shared' if at least two of those pairs have
    no group in common.
    """

    sets = [set(PAIR_RE.match(p).groups()) for p in pairs]

    common = set.intersection(*sets)

    if common:
        return "group-driven"

    for left, right in itertools.combinations(sets, 2):
        if not (left & right):
            return "shared"

    return "group-driven"


def recurrence_null(
    counts_per_pair: dict[str, int],
    total_windows: int,
    replicates: int,
    rng
) -> pd.DataFrame:
    """
    Scatter each pair's outlier windows uniformly at random over the
    genome, repeat, and record how many windows end up flagged by
    exactly k pairs.
    """

    n_pairs = len(counts_per_pair)

    observed_null = np.zeros((replicates, n_pairs + 1), dtype=int)

    for replicate in range(replicates):

        hits = np.zeros(total_windows, dtype=int)

        for how_many in counts_per_pair.values():

            chosen = rng.choice(
                total_windows, size=how_many, replace=False
            )

            hits[chosen] += 1

        for k in range(1, n_pairs + 1):
            observed_null[replicate, k] = int((hits == k).sum())

    return pd.DataFrame({
        "n_pairs": range(1, n_pairs + 1),
        "expected_by_chance": observed_null[:, 1:].mean(axis=0),
        "chance_p95": np.percentile(observed_null[:, 1:], 95, axis=0),
        "chance_max": observed_null[:, 1:].max(axis=0),
    })


def main():

    args = parse_args()

    rng = np.random.default_rng(args.seed)

    regions = load_regions(args.manhattan_dir)

    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f"{regions['pair'].nunique()} comparisons, "
          f"{len(regions)} pair-window records, "
          f"{regions['window'].nunique()} distinct windows\n")


    # ---------------------------------------------------------------
    # CHECK 1
    # ---------------------------------------------------------------

    per_window = (
        regions
        .groupby("window")
        .agg(
            n_pairs=("pair", "nunique"),
            pairs=("pair", lambda s: ",".join(sorted(set(s)))),
            chrom=("chrom", "first"),
            start=("start", "first"),
            max_fst=("max_window_fst", "max"),
            min_snps=("n_snps", "min"),
        )
        .reset_index()
    )

    per_window["kind"] = [
        classify(row.split(","))
        for row in per_window["pairs"]
    ]

    per_window = per_window.sort_values(
        ["n_pairs", "max_fst"], ascending=[False, False]
    )

    observed = (
        per_window["n_pairs"]
        .value_counts()
        .sort_index()
        .rename("observed")
    )

    counts_per_pair = (
        regions
        .groupby("pair")["window"]
        .nunique()
        .to_dict()
    )

    null = recurrence_null(
        counts_per_pair,
        args.total_windows,
        args.replicates,
        rng
    )

    comparison = (
        null
        .set_index("n_pairs")
        .join(observed)
        .fillna(0)
    )

    comparison["observed"] = comparison["observed"].astype(int)

    print("CHECK 1 -- recurrence across the 21 comparisons")
    print("(random scatter of the same number of outliers per pair, "
          f"{args.replicates} replicates)\n")

    print(
        comparison
        .loc[comparison[["observed", "expected_by_chance"]].sum(axis=1) > 0]
        .to_string(
            float_format=lambda v: f"{v:8.2f}",
            columns=[
                "observed", "expected_by_chance", "chance_p95", "chance_max"
            ]
        )
    )

    shared = per_window[
        (per_window["n_pairs"] >= 2) & (per_window["kind"] == "shared")
    ]

    driven = per_window[
        (per_window["n_pairs"] >= 2) & (per_window["kind"] == "group-driven")
    ]

    print(f"\nwindows flagged by 2 or more pairs : "
          f"{int((per_window['n_pairs'] >= 2).sum())}")
    print(f"  of them driven by a single group : {len(driven)}")
    print(f"  of them shared across disjoint pairs : {len(shared)}")

    if len(shared):
        print("\ntop shared windows:")
        print(
            shared
            .head(15)[
                ["window", "n_pairs", "max_fst", "min_snps", "pairs"]
            ]
            .to_string(index=False)
        )

    per_window.to_csv(
        args.outdir / "window_recurrence.tsv", sep="\t", index=False
    )

    comparison.to_csv(
        args.outdir / "recurrence_vs_chance.tsv", sep="\t"
    )


    # ---------------------------------------------------------------
    # CHECK 2
    # ---------------------------------------------------------------

    print("\n\nCHECK 2 -- SNP count against window FST, within outliers\n")

    single = regions[regions["n_windows"] == 1].copy()

    rho, pvalue = stats.spearmanr(
        single["n_snps"], single["mean_window_fst"]
    )

    print(f"single-window outliers        : {len(single)}")
    print(f"SNPs per window, median       : {single['n_snps'].median():.0f}")
    print(f"SNPs per window, 5th pct      : "
          f"{np.percentile(single['n_snps'], 5):.0f}")
    print(f"Spearman(n_snps, window FST)  : {rho:+.3f}  (p = {pvalue:.2g})")

    bins = pd.qcut(single["n_snps"], 4, duplicates="drop")

    by_bin = (
        single
        .groupby(bins, observed=True)
        .agg(
            n=("n_snps", "size"),
            median_snps=("n_snps", "median"),
            mean_fst=("mean_window_fst", "mean"),
            max_fst=("mean_window_fst", "max"),
        )
    )

    print("\nby SNP-count quartile:")
    print(by_bin.to_string(float_format=lambda v: f"{v:8.3f}"))

    single.to_csv(
        args.outdir / "outlier_windows_snp_density.tsv",
        sep="\t",
        index=False
    )

    print(f"\nwritten to {args.outdir}")


if __name__ == "__main__":
    main()
