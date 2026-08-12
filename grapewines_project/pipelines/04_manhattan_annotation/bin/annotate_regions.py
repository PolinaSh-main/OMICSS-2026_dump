#!/usr/bin/env python3
"""
Find the genes that sit in (or next to) each candidate region.

Input:
    <comparison>_candidate_regions.tsv   from merge_regions.py
    PN40024.v4.1.REF.gff3                reference annotation
    PN40024.v4.1.REF.b2g.tsv             Blast2GO functional descriptions

Output:
    <comparison>_candidate_genes.tsv

    region_id  chrom  region_start  region_end  mean_window_fst
    max_window_fst  gene_id  gene_start  gene_end  strand
    distance_to_region_bp  gene_name  description  go_terms

distance_to_region_bp is 0 for a gene that overlaps the region, and the
gap in bp for a gene that only lies within --flank of it. Genes just
outside a sweep are worth keeping -- the causal variant is not
necessarily inside the window that happened to pass the FST threshold --
but they must be distinguishable from the ones actually in it.

Chromosome naming is the usual trap here: the VCF (and therefore the
FST output) uses 1..19, while the reference GFF3 uses chr01..chr19.
Both are reduced to their digits before matching, and the script stops
with a clear error if that still leaves no overlap at all.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path
from urllib.parse import unquote

import pandas as pd


GENE_FEATURES = {"gene", "protein_coding_gene", "ncRNA_gene"}


def parse_args():

    parser = argparse.ArgumentParser(
        description="Annotate candidate FST regions with genes"
    )

    parser.add_argument("--regions",  required=True, type=Path)
    parser.add_argument("--gff",      required=True, type=Path)
    parser.add_argument("--b2g",      default=None,  type=Path,
                        help="Blast2GO table; optional but recommended")
    parser.add_argument("--flank",    type=int, default=10_000,
                        help="Also report genes within this distance of a "
                             "region (default 10000)")
    parser.add_argument("--comparison", default=None)
    parser.add_argument("--outdir",   required=True, type=Path)

    return parser.parse_args()


def normalise_chrom(value) -> str:
    """
    chr01, Chr1, 1, NC_0001.1 -> "1" where possible.

    Anything without digits is lowercased and returned as-is, so
    scaffolds and organelles still match each other.
    """

    text = str(value).strip()

    digits = re.sub(r"\D", "", text)

    if digits:
        return str(int(digits))

    return text.lower()


def open_maybe_gzip(path: Path):

    if path.suffix == ".gz":
        return gzip.open(path, "rt", errors="ignore")

    return path.open("rt", errors="ignore")


def parse_attributes(field: str) -> dict:

    attributes = {}

    for part in field.strip().rstrip(";").split(";"):

        if "=" not in part:
            continue

        key, _, value = part.partition("=")

        attributes[key.strip()] = unquote(value.strip())

    return attributes


def read_genes(path: Path) -> pd.DataFrame:

    rows = []

    with open_maybe_gzip(path) as handle:

        for line in handle:

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                continue

            if fields[2] not in GENE_FEATURES:
                continue

            attributes = parse_attributes(fields[8])

            gene_id = (
                attributes.get("ID")
                or attributes.get("gene_id")
                or attributes.get("Name")
            )

            if gene_id is None:
                continue

            # GFF3 IDs are often prefixed, e.g. "gene:Vitvi01g00001"
            gene_id = gene_id.split(":")[-1]

            rows.append(
                {
                    "chrom_raw": fields[0],
                    "chrom": normalise_chrom(fields[0]),
                    "gene_start": int(fields[3]),
                    "gene_end": int(fields[4]),
                    "strand": fields[6],
                    "gene_id": gene_id,
                    "gene_name": attributes.get("Name", ""),
                    "gff_description": (
                        attributes.get("description")
                        or attributes.get("Note")
                        or ""
                    ),
                }
            )

    if not rows:
        raise ValueError(
            f"No gene features found in {path}. "
            f"Looked for column 3 in {sorted(GENE_FEATURES)}."
        )

    return pd.DataFrame(rows)


def read_b2g(path: Path) -> pd.DataFrame:
    """
    Blast2GO exports vary between versions, so the columns are found by
    name rather than by position. Falls back to "first column is the ID,
    second is the description".
    """

    b2g = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        engine="python",
        on_bad_lines="skip"
    )

    columns = list(b2g.columns)

    def find(*needles, exclude=(), default=None):
        """
        Needles are tried in order of preference, each against every
        column. Scanning columns first instead would let a weak needle
        win: "name" matches the "seq_name" ID column long before it
        reaches a real description column.
        """

        for needle in needles:

            for column in columns:

                if column in exclude:
                    continue

                if needle in column.lower():
                    return column

        return default

    id_column = find(
        "seq", "transcript", "gene", "id",
        default=columns[0]
    )

    description_column = find(
        "description", "product", "annot", "name",
        exclude={id_column},
        default=columns[1] if len(columns) > 1 else None
    )

    go_column = find(
        "go",
        exclude={id_column, description_column},
        default=None
    )

    out = pd.DataFrame(
        {
            "gene_id": b2g[id_column].astype(str).str.split(":").str[-1],
            "description": (
                b2g[description_column]
                if description_column else ""
            ),
            "go_terms": b2g[go_column] if go_column else "",
        }
    )

    # A gene can appear once per transcript; keep the first description.
    return out.drop_duplicates(subset="gene_id")


def strip_transcript_suffix(gene_id: str) -> str:
    """
    Vitvi01g00001_t001 / Vitvi01g00001.1 -> Vitvi01g00001
    """

    return re.sub(r"([._]t?\d+)$", "", str(gene_id))


def main():

    args = parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    comparison = args.comparison or args.regions.name.replace(
        "_candidate_regions.tsv", ""
    )


    regions = pd.read_csv(args.regions, sep="\t")

    if regions.empty:
        print("No candidate regions; nothing to annotate.", file=sys.stderr)

        pd.DataFrame(
            columns=[
                "region_id", "chrom", "region_start", "region_end",
                "mean_window_fst", "max_window_fst", "gene_id",
                "gene_start", "gene_end", "strand",
                "distance_to_region_bp", "gene_name", "description",
                "go_terms",
            ]
        ).to_csv(
            args.outdir / f"{comparison}_candidate_genes.tsv",
            sep="\t",
            index=False
        )

        return

    regions["chrom_key"] = regions["chrom"].map(normalise_chrom)


    genes = read_genes(args.gff)

    print(
        f"{len(genes)} genes in {args.gff.name}, "
        f"{genes['chrom'].nunique()} sequences"
    )


    shared = set(regions["chrom_key"]) & set(genes["chrom"])

    if not shared:
        raise ValueError(
            "No chromosome name is shared between the FST regions and "
            f"the GFF3 after normalisation.\n"
            f"  regions: {sorted(set(regions['chrom_key']))[:10]}\n"
            f"  gff3   : {sorted(set(genes['chrom']))[:10]}"
        )


    hits = []

    for region in regions.itertuples():

        on_chrom = genes[genes["chrom"] == region.chrom_key]

        if on_chrom.empty:
            continue

        # Overlap, or within --flank of either edge.
        near = on_chrom[
            (on_chrom["gene_end"] >= region.start - args.flank)
            & (on_chrom["gene_start"] <= region.end + args.flank)
        ]

        for gene in near.itertuples():

            if gene.gene_end < region.start:
                distance = region.start - gene.gene_end

            elif gene.gene_start > region.end:
                distance = gene.gene_start - region.end

            else:
                distance = 0

            hits.append(
                {
                    "region_id": region.region_id,
                    "chrom": region.chrom,
                    "region_start": region.start,
                    "region_end": region.end,
                    "n_snps": getattr(region, "n_snps", ""),
                    "mean_window_fst": region.mean_window_fst,
                    "max_window_fst": region.max_window_fst,
                    "gene_id": gene.gene_id,
                    "gene_start": gene.gene_start,
                    "gene_end": gene.gene_end,
                    "strand": gene.strand,
                    "distance_to_region_bp": distance,
                    "gene_name": gene.gene_name,
                    "gff_description": gene.gff_description,
                }
            )


    annotated = pd.DataFrame(hits)

    if annotated.empty:
        print(
            "WARNING: candidate regions contain no annotated genes",
            file=sys.stderr
        )


    #
    # Functional descriptions. The GFF3 rarely carries them; Blast2GO
    # is where the readable annotation lives.
    #

    if args.b2g and args.b2g.exists() and not annotated.empty:

        b2g = read_b2g(args.b2g)

        print(f"{len(b2g)} functional annotations in {args.b2g.name}")

        annotated = annotated.merge(b2g, on="gene_id", how="left")

        # Retry the ones that failed, ignoring transcript suffixes.
        unmatched = annotated["description"].isna()

        if unmatched.any():

            b2g_stripped = b2g.assign(
                gene_key=b2g["gene_id"].map(strip_transcript_suffix)
            ).drop_duplicates(subset="gene_key")

            retry = (
                annotated.loc[unmatched, ["gene_id"]]
                .assign(
                    gene_key=lambda d: d["gene_id"].map(
                        strip_transcript_suffix
                    )
                )
                .merge(
                    b2g_stripped[["gene_key", "description", "go_terms"]],
                    on="gene_key",
                    how="left"
                )
            )

            annotated.loc[unmatched, "description"] = (
                retry["description"].values
            )

            annotated.loc[unmatched, "go_terms"] = (
                retry["go_terms"].values
            )

        print(
            f"{annotated['description'].notna().sum()} of "
            f"{len(annotated)} gene hits matched a description"
        )

    else:

        annotated["description"] = annotated.get("gff_description", "")
        annotated["go_terms"] = ""


    # Fall back to whatever the GFF3 said when Blast2GO had nothing.
    annotated["description"] = (
        annotated["description"]
        .fillna("")
        .replace("", pd.NA)
        .fillna(annotated["gff_description"])
        .fillna("")
    )

    annotated = annotated.drop(columns=["gff_description"])


    annotated = annotated.sort_values(
        ["max_window_fst", "mean_window_fst", "distance_to_region_bp"],
        ascending=[False, False, True]
    ).reset_index(drop=True)


    output = args.outdir / f"{comparison}_candidate_genes.tsv"

    annotated.to_csv(
        output,
        sep="\t",
        index=False,
        float_format="%.6f"
    )

    print()
    print(
        f"{len(annotated)} gene hits across "
        f"{annotated['region_id'].nunique() if not annotated.empty else 0} "
        f"regions -> {output}"
    )


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
