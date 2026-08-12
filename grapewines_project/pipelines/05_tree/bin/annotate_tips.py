#!/usr/bin/env python3
"""
Annotate the tips of the rooted tree with the K = 7 ADMIXTURE result.

Assignment rule, as set by the task: a sample joins the group where its
ancestry proportion is largest, but only if that proportion reaches
--min-q. Below the threshold it is Admixed. The outgroup is not part of
the K = 7 classification at all and gets its own category.

The group *labels* are derived here from the metadata rather than copied
from any table of K numbers. ADMIXTURE numbers its columns arbitrarily,
so a given K index means whatever the run made it mean; a label taken
from a different run would put confident, wrong text in the legend.

Output:
    tip_annotation.tsv   one row per tip: category, K group, max Q
    group_labels.tsv     what each K group turned out to be
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import newick


ID_COLUMN_CANDIDATES = ["Column 1", "IID", "ID", "id", "sample", "Sample"]


def parse_args():

    parser = argparse.ArgumentParser(
        description="Build the tip annotation table for the rooted tree"
    )

    parser.add_argument("--tree", required=True, type=Path)

    parser.add_argument("--q", required=True, type=Path)

    parser.add_argument(
        "--fam",
        required=True,
        type=Path,
        help="The .fam of the same ADMIXTURE run; supplies the row order"
    )

    parser.add_argument("--metadata", required=True, type=Path)

    parser.add_argument("--min-q", type=float, default=0.75)

    parser.add_argument("--outgroup", default="ZZ01")

    parser.add_argument("--outdir", required=True, type=Path)

    return parser.parse_args()


def normalise(name: str) -> str:
    """
    PLINK writes FID_IID. Collapse only the exact X_X form -- some real
    accession names contain underscores, so a blanket strip would merge
    distinct samples.
    """

    first, sep, second = name.partition("_")

    if sep and first == second:
        return first

    return name


def find_id_column(table: pd.DataFrame) -> str:

    for column in ID_COLUMN_CANDIDATES:

        if column in table.columns:
            return column

    raise ValueError(
        f"No sample ID column in the metadata; saw {list(table.columns)}"
    )


def describe(members: pd.DataFrame) -> tuple[str, str, str]:
    """
    (dominant country, wild or cultivated, a readable label).
    """

    if members.empty:
        return "", "", ""

    countries = members["Geographic origin by country"].dropna()

    country = (
        countries.value_counts().idxmax().title()
        if not countries.empty else "unknown"
    )

    share = (
        countries.value_counts(normalize=True).iloc[0]
        if not countries.empty else 0.0
    )

    background = members["Genetic background"].fillna("").str.lower()

    wild = background.str.contains("sylvestris").sum()

    cultivated = background.str.contains("ssp. vinifera").sum()

    kind = "wild" if wild > cultivated else "cultivated"

    #
    # A group drawn from two countries in comparable numbers is named
    # after both; anything below two thirds is not "the Armenian group".
    #
    if share < 0.66 and len(countries.value_counts()) > 1:

        second = countries.value_counts().index[1].title()

        country = f"{country}-{second}"

    label = f"{country} {kind}"

    return country, kind, label


def main():

    args = parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)


    tree = newick.parse(args.tree.read_text())

    tips = [normalise(name) for name in newick.tip_names(tree)]

    if len(set(tips)) != len(tips):

        seen, duplicated = set(), sorted(
            {name for name in tips if tips.count(name) > 1}
        )

        raise ValueError(
            f"Duplicated tip labels in the tree: {', '.join(duplicated[:10])}"
        )

    print(f"{len(tips)} tips")


    #
    # Q rows are in .fam order. The two files must come from the same
    # ADMIXTURE run; a mismatched pair is silently wrong, so the row
    # counts are checked rather than trusted.
    #

    fam = pd.read_csv(
        args.fam,
        sep=r"\s+",
        header=None,
        dtype=str,
        usecols=[0, 1],
        names=["FID", "IID"],
    )

    q = pd.read_csv(args.q, sep=r"\s+", header=None)

    if len(fam) != len(q):
        raise ValueError(
            f"{args.fam.name} has {len(fam)} rows but {args.q.name} has "
            f"{len(q)}; they are not from the same ADMIXTURE run"
        )

    n_components = q.shape[1]

    print(f"{len(q)} samples x K={n_components}")

    sample_ids = [normalise(name) for name in fam["IID"]]

    values = q.to_numpy(dtype=float)

    max_q = values.max(axis=1)

    best = values.argmax(axis=1) + 1

    assignment = pd.DataFrame(
        {
            "sample": sample_ids,
            "k_group": [f"K{index}" for index in best],
            "max_q": max_q,
        }
    )

    assignment["category"] = np.where(
        assignment["max_q"] >= args.min_q,
        assignment["k_group"],
        "Admixed",
    )


    #
    # Every tip needs exactly one annotation, and every annotated sample
    # needs to be in the tree.
    #

    by_sample = assignment.set_index("sample")

    missing = [
        tip for tip in tips
        if tip not in by_sample.index and tip != args.outgroup
    ]

    extra = [
        sample for sample in assignment["sample"]
        if sample not in set(tips)
    ]

    if missing or extra:

        raise ValueError(
            "Tree and annotation do not line up.\n"
            f"  in the tree but not annotated ({len(missing)}): "
            f"{', '.join(missing[:10])}\n"
            f"  annotated but not in the tree ({len(extra)}): "
            f"{', '.join(extra[:10])}"
        )


    rows = []

    for tip in tips:

        if tip == args.outgroup:

            rows.append(
                {
                    "tip": tip,
                    "category": "Outgroup",
                    "k_group": "",
                    "max_q": np.nan,
                }
            )

            continue

        record = by_sample.loc[tip]

        rows.append(
            {
                "tip": tip,
                "category": record["category"],
                "k_group": record["k_group"],
                "max_q": round(float(record["max_q"]), 4),
            }
        )

    annotation = pd.DataFrame(rows)

    annotation.to_csv(
        args.outdir / "tip_annotation.tsv",
        sep="\t",
        index=False,
    )


    #
    # What each group actually is, from the metadata.
    #

    metadata = pd.read_csv(args.metadata)

    id_column = find_id_column(metadata)

    metadata[id_column] = metadata[id_column].astype(str).map(normalise)

    metadata = metadata.set_index(id_column)

    label_rows = []

    for group in [f"K{index}" for index in range(1, n_components + 1)]:

        members = annotation.loc[annotation["category"] == group, "tip"]

        present = [tip for tip in members if tip in metadata.index]

        country, kind, label = describe(metadata.loc[present])

        label_rows.append(
            {
                "group": group,
                "n": len(members),
                "main_country": country,
                "wild_or_cultivated": kind,
                "label": label,
            }
        )

    #
    # Two K groups can describe to the same words -- there really are
    # two distinct Armenian wild groups here, and a legend with the same
    # text twice is useless. Number them in K order.
    #

    seen: dict[str, int] = {}

    for row in label_rows:

        if row["label"]:
            seen[row["label"]] = seen.get(row["label"], 0) + 1

    running: dict[str, int] = {}

    for row in label_rows:

        if row["label"] and seen[row["label"]] > 1:

            running[row["label"]] = running.get(row["label"], 0) + 1

            row["label"] = f"{row['label']} group {running[row['label']]}"

    label_rows.append(
        {
            "group": "Admixed",
            "n": int((annotation["category"] == "Admixed").sum()),
            "main_country": "",
            "wild_or_cultivated": "",
            "label": f"Admixed (max Q < {args.min_q})",
        }
    )

    label_rows.append(
        {
            "group": "Outgroup",
            "n": int((annotation["category"] == "Outgroup").sum()),
            "main_country": "",
            "wild_or_cultivated": "",
            "label": f"{args.outgroup} (V. rotundifolia)",
        }
    )

    labels = pd.DataFrame(label_rows)

    labels.to_csv(
        args.outdir / "group_labels.tsv",
        sep="\t",
        index=False,
    )

    print()
    print(labels.to_string(index=False))
    print()
    print(f"{int((annotation['category'] == 'Admixed').sum())} tips admixed, "
          f"{int((annotation['category'] != 'Admixed').sum()) - 1} assigned")


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
