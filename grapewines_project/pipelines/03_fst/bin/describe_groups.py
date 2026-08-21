#!/usr/bin/env python3
"""
Turn anonymous ADMIXTURE groups into named, interpretable ones.

ADMIXTURE only says "component 5"; the report needs "Armenian wild
group". This script joins the per-sample assignment table to the
metadata and works out, for every group, which country and which
domestication status dominate it.

Input:
    sample_q_values_K<K>.tsv   from assign_admixture_groups.py
    metadata CSV               one row per sample

Output:
    group_composition_K<K>.tsv    every group x every metadata level
    group_interpretation_K<K>.tsv one row per group: proposed name + evidence
    group_labels_K<K>.tsv         group -> short label, for figure axes

Nothing here is authoritative -- it is a first pass over the metadata
that still has to be read and confirmed before it goes into a report.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


#
# Metadata columns worth breaking every group down by.
#

PROFILE_COLUMNS = [
    "Geographic origin by country",
    "Geographic origin by region",
    "Genetic background",
    "Utilization",
    "Berry skin color",
]


ID_COLUMN_CANDIDATES = [
    "Column 1",
    "IID",
    "ID",
    "id",
    "sample",
    "Sample",
]


#
# A group is named after a single country only when effectively all of
# its samples come from there; otherwise it is "<COUNTRY>-dominated ...
# (heterogeneous)". At K=7 this puts K1, K2, K5 and K7 (all at 100%) in
# the first form and K3, K4 and K6 (46%, 58%, 73%) in the second, which
# is the wording used on the group figures.
#
HOMOGENEOUS_FRACTION = 0.95


def parse_args():

    parser = argparse.ArgumentParser(
        description="Describe and name ADMIXTURE groups from metadata"
    )

    parser.add_argument(
        "--assignments",
        required=True,
        type=Path,
        help="sample_q_values_K<K>.tsv"
    )

    parser.add_argument(
        "--metadata",
        required=True,
        type=Path
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path
    )

    return parser.parse_args()


def find_id_column(metadata: pd.DataFrame) -> str:

    for column in ID_COLUMN_CANDIDATES:

        if column in metadata.columns:
            return column

    raise ValueError(
        f"No sample ID column in metadata. "
        f"Available: {list(metadata.columns)}"
    )


def classify_domestication(row) -> str:
    """
    wild / cultivated / other, from whichever column answers first.

    ssp. sylvestris is the wild grapevine, ssp. vinifera the cultivated
    one; Utilization is the fallback for samples whose subspecies is
    blank or a cross.
    """

    background = str(row.get("Genetic background", "")).lower()

    if "sylvestris" in background:
        return "wild"

    if "vinifera ssp. vinifera" in background:
        return "cultivated"

    utilization = str(row.get("Utilization", "")).lower()

    if utilization in ("wild", "feral"):
        return "wild"

    if utilization in (
        "wine", "table", "raisin",
        "wine/table", "table/raisin", "wine/table/raisin"
    ):
        return "cultivated"

    return "other"


def top_levels(series: pd.Series, drop_na: bool = True):
    """
    Value counts as (level, n, fraction), most common first.
    """

    values = series.astype(str).str.strip()

    if drop_na:
        values = values[~values.isin(["NA", "nan", ""])]

    if len(values) == 0:
        return []

    counts = values.value_counts()

    return [
        (level, int(n), n / len(values))
        for level, n in counts.items()
    ]


def propose_name(country_levels, domestication_levels) -> tuple[str, str]:
    """
    Returns (name, evidence).

    The wording follows the naming convention agreed for the group
    figures:

        ARMENIA wild
        TURKEY-dominated cultivated (heterogeneous)

    A group is named after one country only when effectively all of its
    samples come from there. Below that it is "<COUNTRY>-dominated ...
    (heterogeneous)", which says both where the group sits and that the
    label is a summary rather than a fact about every member.

    Two groups can end up with the same words -- there really are two
    all-Armenian wild groups here. They are not numbered apart: the K
    index is already in front of the label everywhere it is used, and
    inventing "group 1"/"group 2" would imply an ordering the data does
    not support.
    """

    if not country_levels:
        return "unassigned group", "no country metadata"

    country, n_country, frac_country = country_levels[0]

    if domestication_levels:

        status, _, frac_status = domestication_levels[0]

        status_evidence = f"{frac_status:.0%} {status}"

    else:

        status = "unclassified"

        status_evidence = "no domestication metadata"

    if frac_country >= HOMOGENEOUS_FRACTION:

        name = f"{country.upper()} {status}"

        country_evidence = (
            f"{frac_country:.0%} of samples from {country}"
        )

    else:

        second = (
            country_levels[1][0]
            if len(country_levels) > 1
            else None
        )

        name = f"{country.upper()}-dominated {status} (heterogeneous)"

        country_evidence = (
            f"mixed origin, top countries {country} "
            f"({frac_country:.0%})"
            + (
                f" and {second} ({country_levels[1][2]:.0%})"
                if second else ""
            )
        )

    return name, f"{country_evidence}; {status_evidence}"


def main():

    args = parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)


    assignments = pd.read_csv(
        args.assignments,
        sep="\t",
        dtype={"FID": str, "IID": str}
    )

    metadata = pd.read_csv(args.metadata)

    id_column = find_id_column(metadata)

    metadata[id_column] = metadata[id_column].astype(str)


    merged = assignments.merge(
        metadata,
        left_on="IID",
        right_on=id_column,
        how="left"
    )


    matched = merged[id_column].notna().sum()

    if matched == 0:
        raise ValueError(
            "No sample ID in the .fam matches the metadata"
        )

    if matched < len(merged):
        print(
            f"WARNING: {len(merged) - matched} of {len(merged)} samples "
            f"have no metadata row",
            file=sys.stderr
        )


    merged["domestication"] = merged.apply(
        classify_domestication,
        axis=1
    )


    #
    # K is the number of ancestry columns, named K1..Kn upstream.
    #

    group_columns = [
        c for c in assignments.columns
        if re.fullmatch(r"K\d+", c)
    ]

    k = len(group_columns)


    present = [
        g for g in group_columns
        if (merged["assignment"] == g).any()
    ]


    composition_rows = []

    interpretation_rows = []


    for group in present:

        members = merged[merged["assignment"] == group]


        for column in PROFILE_COLUMNS + ["domestication"]:

            if column not in members.columns:
                continue

            for level, n, fraction in top_levels(members[column]):

                composition_rows.append(
                    {
                        "group": group,
                        "variable": column,
                        "level": level,
                        "n_samples": n,
                        "fraction": round(fraction, 4),
                    }
                )


        countries = top_levels(
            members["Geographic origin by country"]
        )

        domestication = top_levels(
            members.loc[
                members["domestication"] != "other",
                "domestication"
            ]
        )

        name, evidence = propose_name(countries, domestication)


        interpretation_rows.append(
            {
                "group": group,
                "n_samples": len(members),
                "main_country": countries[0][0] if countries else "NA",
                "main_country_fraction": (
                    round(countries[0][2], 3) if countries else float("nan")
                ),
                "second_country": (
                    countries[1][0] if len(countries) > 1 else ""
                ),
                "wild_or_cultivated": (
                    domestication[0][0] if domestication else "NA"
                ),
                "domestication_fraction": (
                    round(domestication[0][2], 3)
                    if domestication else float("nan")
                ),
                "homogeneous": (
                    "yes"
                    if countries and countries[0][2] >= HOMOGENEOUS_FRACTION
                    else "no"
                ),
                "proposed_name": name,
                "evidence": evidence,
            }
        )


    composition = pd.DataFrame(composition_rows)

    interpretation = pd.DataFrame(interpretation_rows)


    composition.to_csv(
        args.outdir / f"group_composition_K{k}.tsv",
        sep="\t",
        index=False
    )

    interpretation.to_csv(
        args.outdir / f"group_interpretation_K{k}.tsv",
        sep="\t",
        index=False
    )


    #
    # Short axis labels for the heatmap: "K5 - Armenian wild group 1"
    #

    labels = interpretation[["group", "proposed_name"]].copy()

    labels["label"] = (
        labels["group"] + " - " + labels["proposed_name"]
    )

    labels[["group", "label"]].to_csv(
        args.outdir / f"group_labels_K{k}.tsv",
        sep="\t",
        index=False
    )


    print(
        interpretation[
            [
                "group",
                "n_samples",
                "main_country",
                "wild_or_cultivated",
                "proposed_name",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
