#!/usr/bin/env python3
"""
Collect the equal-n FST replicates and compare them with the full-group
run from 03_fst.

Answers three questions:

  1. Does the size/FST correlation survive equalising the sample sizes?
     Spearman between the ORIGINAL group sizes and the subsampled FST.
     If it collapses towards zero, the ranking in 03_fst was a sampling
     artefact. If it holds, the narrow groups are genuinely narrow.

  2. How much does the estimate move between draws? A pair whose FST
     wanders by more than the gaps between pairs cannot be ranked at
     all, whatever the mean says.

  3. Do the groups differ in within-group diversity? pi is literally the
     denominator of FST, so this is the direct measurement rather than
     the inference.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


MEAN_RE = re.compile(
    r"Weir and Cockerham mean Fst estimate:\s*([0-9.eE+-]+)"
)

WEIGHTED_RE = re.compile(
    r"Weir and Cockerham weighted Fst estimate:\s*([0-9.eE+-]+)"
)

FST_NAME_RE = re.compile(r"^rep(\d+)_(K\d+)_vs_(K\d+)$")

PI_NAME_RE = re.compile(r"^(rep\d+|full)_(K\d+)\.pi$")


def parse_args():

    parser = argparse.ArgumentParser(
        description="Summarise equal-n FST replicates"
    )

    parser.add_argument("--fstdir", required=True, type=Path)

    parser.add_argument("--pidir", required=True, type=Path)

    parser.add_argument(
        "--full-summary",
        required=True,
        type=Path,
        help="pairwise_fst_summary_K<K>.tsv from 03_fst"
    )

    parser.add_argument(
        "--group-sizes",
        required=True,
        type=Path,
        help="group_sizes_K<K>.tsv from 03_fst"
    )

    parser.add_argument("--outdir", required=True, type=Path)

    return parser.parse_args()


def read_replicates(fstdir: Path) -> pd.DataFrame:

    rows = []

    for path in sorted(fstdir.glob("rep*_K*_vs_K*.log")):

        match = FST_NAME_RE.match(path.stem)

        if match is None:
            continue

        replicate, group_a, group_b = match.groups()

        text = path.read_text(errors="ignore")

        mean = MEAN_RE.search(text)

        weighted = WEIGHTED_RE.search(text)

        if mean is None:
            raise ValueError(f"No mean FST in {path}")

        rows.append({
            "replicate": int(replicate),
            "comparison": f"{group_a}_vs_{group_b}",
            "group_a": group_a,
            "group_b": group_b,
            "mean_fst": float(mean.group(1)),
            "weighted_fst": (
                float(weighted.group(1)) if weighted else np.nan
            ),
        })

    if not rows:
        raise FileNotFoundError(f"No rep*_K*_vs_K*.log in {fstdir}")

    return pd.DataFrame(rows)


def read_pi(pidir: Path) -> pd.DataFrame:

    rows = []

    for path in sorted(pidir.glob("*_K*.pi")):

        match = PI_NAME_RE.match(path.name)

        if match is None:
            continue

        scope, group = match.groups()

        value = path.read_text().strip()

        rows.append({
            "scope": "full" if scope == "full" else "equal_n",
            "replicate": 0 if scope == "full" else int(scope[3:]),
            "group": group,
            "pi": float(value) if value else np.nan,
        })

    return pd.DataFrame(rows)


def main():

    args = parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    replicates = read_replicates(args.fstdir)

    full = pd.read_csv(args.full_summary, sep="\t")

    sizes = pd.read_csv(args.group_sizes, sep="\t")

    sizes = dict(
        zip(sizes["group"], sizes["n_samples"])
    )

    n_replicates = replicates["replicate"].nunique()


    # ---------------------------------------------------------------
    # per-comparison summary
    # ---------------------------------------------------------------

    summary = (
        replicates
        .groupby(["comparison", "group_a", "group_b"])
        .agg(
            n_replicates=("mean_fst", "size"),
            equal_n_mean=("mean_fst", "mean"),
            equal_n_sd=("mean_fst", "std"),
            equal_n_min=("mean_fst", "min"),
            equal_n_max=("mean_fst", "max"),
        )
        .reset_index()
    )

    summary = summary.merge(
        full[["comparison", "mean_fst"]].rename(
            columns={"mean_fst": "full_group_fst"}
        ),
        on="comparison",
        how="left"
    )

    summary["size_a"] = summary["group_a"].map(sizes)

    summary["size_b"] = summary["group_b"].map(sizes)

    summary["smaller_group"] = summary[["size_a", "size_b"]].min(axis=1)

    summary["mean_group_size"] = (
        summary[["size_a", "size_b"]].mean(axis=1)
    )

    summary["change"] = (
        summary["equal_n_mean"] - summary["full_group_fst"]
    )

    summary = summary.sort_values("equal_n_mean", ascending=False)


    # ---------------------------------------------------------------
    # question 1 -- does the size effect survive?
    # ---------------------------------------------------------------

    def spearman(column, against):
        rho, p = stats.spearmanr(summary[against], summary[column])
        return rho, p

    before_rho, before_p = spearman("full_group_fst", "mean_group_size")

    after_rho, after_p = spearman("equal_n_mean", "mean_group_size")

    before_small, _ = spearman("full_group_fst", "smaller_group")

    after_small, _ = spearman("equal_n_mean", "smaller_group")

    print("=" * 66)
    print("EQUAL-n FST  --  {} replicates, {} comparisons".format(
        n_replicates, summary["comparison"].nunique()
    ))
    print("=" * 66)

    print("\nQ1  Spearman(group size, FST)\n")
    print(f"    full groups, mean size    : {before_rho:+.3f}  "
          f"(p = {before_p:.3g})")
    print(f"    equal n,     mean size    : {after_rho:+.3f}  "
          f"(p = {after_p:.3g})")
    print(f"    full groups, smaller of 2 : {before_small:+.3f}")
    print(f"    equal n,     smaller of 2 : {after_small:+.3f}")


    # ---------------------------------------------------------------
    # question 2 -- is the ranking stable?
    # ---------------------------------------------------------------

    spread = summary["equal_n_sd"].median()

    gaps = np.diff(np.sort(summary["equal_n_mean"].values))

    print("\nQ2  is the ranking resolvable?\n")
    print(f"    median SD between draws   : {spread:.4f}")
    print(f"    median gap between pairs  : {np.median(gaps):.4f}")
    print(f"    -> {'ranking is finer than the noise; do not rank' if spread > np.median(gaps) else 'gaps exceed the noise; ranking is readable'}")

    rank_before = summary["full_group_fst"].rank(ascending=False)

    rank_after = summary["equal_n_mean"].rank(ascending=False)

    tau, tau_p = stats.kendalltau(rank_before, rank_after)

    print(f"    Kendall tau, old vs new   : {tau:+.3f} (p = {tau_p:.3g})")

    print("\n" + summary[[
        "comparison", "size_a", "size_b",
        "full_group_fst", "equal_n_mean", "equal_n_sd", "change"
    ]].to_string(index=False, float_format=lambda v: f"{v:9.4f}"))


    # ---------------------------------------------------------------
    # question 3 -- within-group diversity
    # ---------------------------------------------------------------

    pi = read_pi(args.pidir)

    if not pi.empty:

        print("\n\nQ3  within-group nucleotide diversity (the denominator)\n")

        wide = (
            pi
            .groupby(["group", "scope"])["pi"]
            .agg(["mean", "std"])
            .reset_index()
            .pivot(index="group", columns="scope", values=["mean", "std"])
        )

        wide.columns = [f"{b}_{a}" for a, b in wide.columns]

        wide["n_samples"] = wide.index.map(sizes)

        wide = wide.sort_values("n_samples")

        print(wide.to_string(float_format=lambda v: f"{v:10.5f}"))

        full_pi = pi[pi["scope"] == "full"]

        if len(full_pi) > 2:

            rho, p = stats.spearmanr(
                full_pi["group"].map(sizes), full_pi["pi"]
            )

            print(f"\n    Spearman(group size, pi) on full groups : "
                  f"{rho:+.3f} (p = {p:.3g})")

            print("    Positive means the larger groups carry more "
                  "diversity, which is\n    what links size to FST. "
                  "With seven groups this has almost no\n    power on "
                  "its own -- Q4 below is the test that decides it.")

        pi.to_csv(args.outdir / "within_group_pi.tsv", sep="\t", index=False)

        # -----------------------------------------------------------
        # question 4 -- is FST here just a restatement of Hs?
        # -----------------------------------------------------------
        #
        # FST = 1 - Hs/Ht, where Hs is diversity within the groups and
        # Ht diversity in the two of them pooled. If Ht turns out to be
        # the same for every pair, then nothing in the ranking is about
        # how far apart two groups are -- the whole spread comes from
        # Hs, i.e. from how uniform the groups are internally.
        #
        # Worth checking because a 12% spread in Hs does not look like
        # it could produce a threefold spread in FST, yet it can: FST is
        # one minus a ratio close to one, so small moves in Hs are
        # amplified.
        #

        pi_full = dict(
            zip(full_pi["group"], full_pi["pi"])
        )

        summary["Hs"] = (
            summary["group_a"].map(pi_full)
            + summary["group_b"].map(pi_full)
        ) / 2

        summary["implied_Ht"] = (
            summary["Hs"] / (1 - summary["equal_n_mean"])
        )

        hs_rho, hs_p = stats.spearmanr(
            summary["Hs"], summary["equal_n_mean"]
        )

        print("\n\nQ4  is the ranking a restatement of within-group "
              "diversity?\n")

        print(f"    Spearman(Hs of the pair, FST) : {hs_rho:+.3f} "
              f"(p = {hs_p:.2g})")

        print(f"    Spearman(size of pair,   FST) : {after_rho:+.3f}")

        print(f"\n    implied Ht = Hs / (1 - FST)   : "
              f"{summary['implied_Ht'].mean():.4f} "
              f"+/- {summary['implied_Ht'].std():.4f}   "
              f"({summary['implied_Ht'].min():.4f} to "
              f"{summary['implied_Ht'].max():.4f})")

        constant_ht = summary["implied_Ht"].mean()

        print("\n    holding Ht at that mean and moving only Hs:")

        for label, value in [
            ("least diverse pair", summary["Hs"].min()),
            ("most diverse pair", summary["Hs"].max()),
        ]:
            print(f"        {label:20s} Hs = {value:.4f} "
                  f"-> FST = {1 - value / constant_ht:.4f}")

        print(f"        {'observed':20s}         "
              f"      -> FST = {summary['equal_n_mean'].min():.4f} "
              f"to {summary['equal_n_mean'].max():.4f}")

    summary.to_csv(
        args.outdir / "fst_equal_n_summary.tsv",
        sep="\t",
        index=False,
        float_format="%.6f"
    )

    replicates.to_csv(
        args.outdir / "fst_equal_n_replicates.tsv",
        sep="\t",
        index=False,
        float_format="%.6f"
    )

    print(f"\nwritten to {args.outdir}")


if __name__ == "__main__":
    main()
