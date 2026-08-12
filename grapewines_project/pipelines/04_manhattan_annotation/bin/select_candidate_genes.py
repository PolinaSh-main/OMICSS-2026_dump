#!/usr/bin/env python3
"""
Pick the genes worth putting in the report.

Taking the first rows of the annotation table does not work: several
genes usually sit inside the same high-FST region, so the top of the
table is one locus repeated. The selection follows the rules agreed for
this project:

  1. rank distinct regions by max_window_fst, then mean_window_fst;
  2. keep the top --top-regions of them;
  3. inside each region prefer genes that actually overlap it
     (distance_to_region_bp = 0);
  4. among those prefer a clear functional annotation -- stress
     response, disease resistance, development, hormone signalling,
     berry metabolism -- and keep one or two;
  5. put hypothetical, unnamed and transposon-related genes in a
     separate list, and fall back to them only when a region has
     nothing better.

Input:
    <comparison>_candidate_genes.tsv   from annotate_regions.py

Output:
    <comparison>_selected_regions.tsv       the report table
    <comparison>_low_information_genes.tsv  the ones set aside
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


#
# Themes the project cares about, in the order they were given. The
# first theme a description matches supplies the note.
#

THEMES = [
    (
        "disease resistance",
        r"\b(resistan\w*|disease|pathogen|defen[cs]e|"
        r"nbs-lrr|nb-arc|lrr|rpp\d|powdery mildew|downy mildew|"
        r"chitinase|thaumatin|pr-?\d)\b",
    ),
    (
        "stress response",
        r"\b(stress|drought|cold|heat shock|hsp\d*|salt|osmotic|"
        r"dehydrin|late embryogenesis|oxidative|wound\w*)\b",
    ),
    (
        "hormone signalling",
        r"\b(auxin|gibberellin|abscisic|aba |ethylene|cytokinin|"
        r"jasmon\w*|salicyl\w*|brassinosteroid|hormone|"
        r"receptor kinase)\b",
    ),
    (
        "berry metabolism",
        r"\b(anthocyanin|flavonoid|stilbene|resveratrol|terpene|"
        r"terpenoid|myb\d*|chalcone|udp-glyco\w*|sugar transport\w*|"
        r"invertase|malate|tartrate|aroma|berry)\b",
    ),
    (
        "development",
        r"\b(development\w*|meristem|floral|flower\w*|"
        r"transcription factor|homeobox|mads|growth|"
        r"cell wall|embryo\w*|seed)\b",
    ),
]


#
# Annotations that carry no biological information.
#

LOW_INFORMATION = re.compile(
    r"("
    r"hypothetical|uncharacteri[sz]ed|unknown|predicted protein|"
    r"putative uncharacteri[sz]ed|"
    r"transpos\w*|retrotranspos\w*|gypsy|copia|helitron|"
    r"integrase|reverse transcriptase|"
    r"^\s*$"
    r")",
    re.IGNORECASE,
)


def parse_args():

    parser = argparse.ArgumentParser(
        description="Select report-worthy genes from annotated regions"
    )

    parser.add_argument("--genes", required=True, type=Path)

    parser.add_argument(
        "--top-regions",
        type=int,
        default=5,
        help="How many distinct regions to report (default 5)"
    )

    parser.add_argument(
        "--genes-per-region",
        type=int,
        default=2,
        help="At most this many genes per region (default 2)"
    )

    parser.add_argument("--comparison", default=None)

    parser.add_argument("--outdir", required=True, type=Path)

    return parser.parse_args()


def classify(description: str, gene_name: str) -> tuple[str, bool]:
    """
    Returns (theme, is_low_information).

    theme is "" when the annotation is readable but does not fall into
    one of the project's themes.
    """

    text = f"{description} {gene_name}".strip()

    if not text or LOW_INFORMATION.search(text):
        return "", True

    for theme, pattern in THEMES:

        if re.search(pattern, text, flags=re.IGNORECASE):
            return theme, False

    return "", False


def rank_regions(genes: pd.DataFrame) -> list[str]:

    regions = (
        genes[["region_id", "max_window_fst", "mean_window_fst"]]
        .drop_duplicates(subset="region_id")
        .sort_values(
            ["max_window_fst", "mean_window_fst"],
            ascending=False
        )
    )

    return regions["region_id"].tolist()


def choose_for_region(
    block: pd.DataFrame,
    limit: int
) -> pd.DataFrame:
    """
    Overlapping first, then themed, then any readable annotation, then
    -- only if there is nothing else -- a low-information gene so the
    region is still represented.
    """

    overlapping = block[block["distance_to_region_bp"] == 0]

    pool = overlapping if not overlapping.empty else block

    themed = pool[
        (pool["theme"] != "") & (~pool["low_information"])
    ]

    readable = pool[~pool["low_information"]]

    for candidates in (themed, readable, pool):

        if not candidates.empty:

            return candidates.sort_values(
                ["distance_to_region_bp", "max_window_fst"],
                ascending=[True, False]
            ).head(limit)

    return block.head(0)


def make_note(row) -> str:

    parts = []

    if row["theme"]:
        parts.append(f"annotation points to {row['theme']}")

    if row["distance_to_region_bp"] == 0:
        parts.append("overlaps the high-FST region")
    else:
        parts.append(
            f"{int(row['distance_to_region_bp'])} bp from the region"
        )

    if row["low_information"]:
        parts.append(
            "no informative annotation in this region -- treat as a "
            "placeholder"
        )

    return "; ".join(parts)


def main():

    args = parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    comparison = args.comparison or args.genes.name.replace(
        "_candidate_genes.tsv", ""
    )


    genes = pd.read_csv(args.genes, sep="\t")

    if genes.empty:
        print("No annotated genes to select from.", file=sys.stderr)

        for suffix in ("selected_regions", "low_information_genes"):
            (args.outdir / f"{comparison}_{suffix}.tsv").write_text(
                "region_id\n"
            )

        return


    for column in ("description", "gene_name"):
        if column not in genes.columns:
            genes[column] = ""

    genes["description"] = genes["description"].fillna("").astype(str)

    genes["gene_name"] = genes["gene_name"].fillna("").astype(str)


    classified = genes.apply(
        lambda row: classify(row["description"], row["gene_name"]),
        axis=1
    )

    genes["theme"] = [theme for theme, _ in classified]

    genes["low_information"] = [low for _, low in classified]


    #
    # Everything uninformative, kept aside as instructed.
    #

    low_information = genes[genes["low_information"]].copy()

    low_information.drop(
        columns=["theme", "low_information"]
    ).to_csv(
        args.outdir / f"{comparison}_low_information_genes.tsv",
        sep="\t",
        index=False,
        float_format="%.6f"
    )


    ranked = rank_regions(genes)[: args.top_regions]

    selected = []

    for rank, region_id in enumerate(ranked, start=1):

        block = genes[genes["region_id"] == region_id]

        chosen = choose_for_region(block, args.genes_per_region)

        for row in chosen.to_dict("records"):
            row["region_rank"] = rank
            selected.append(row)


    if not selected:
        raise RuntimeError(
            "No gene could be selected from the top regions"
        )

    report = pd.DataFrame(selected)

    report["note"] = report.apply(make_note, axis=1)

    columns = [
        "region_rank",
        "region_id",
        "chrom",
        "region_start",
        "region_end",
        "n_snps",
        "mean_window_fst",
        "max_window_fst",
        "gene_id",
        "gene_name",
        "distance_to_region_bp",
        "description",
        "theme",
        "note",
    ]

    report = report[[c for c in columns if c in report.columns]]

    output = args.outdir / f"{comparison}_selected_regions.tsv"

    report.to_csv(
        output,
        sep="\t",
        index=False,
        float_format="%.6f"
    )


    print(
        f"{len(ranked)} regions reported, "
        f"{len(report)} genes selected, "
        f"{len(low_information)} genes set aside"
    )

    print()

    print(
        report[
            [
                "region_rank", "region_id", "gene_id",
                "max_window_fst", "theme"
            ]
        ].to_string(index=False)
    )

    print()
    print(output)


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
