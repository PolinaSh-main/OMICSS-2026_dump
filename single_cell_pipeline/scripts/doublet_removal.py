#!/usr/bin/env python3

import argparse
from pathlib import Path

import scanpy as sc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)

    parser.add_argument(
        "--expected-doublet-rate",
        type=float,
        default=0.06,
    )

    args = parser.parse_args()

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.input)

    # Scrublet expects raw counts.
    sc.pp.scrublet(
        adata,
        expected_doublet_rate=args.expected_doublet_rate,
    )

    # Save all predictions before filtering
    adata.obs[
        [
            "doublet_score",
            "predicted_doublet",
        ]
    ].to_csv(
        output / "doublet_predictions.tsv",
        sep="\t",
    )

    n_before = adata.n_obs

    adata = adata[
        ~adata.obs["predicted_doublet"]
    ].copy()

    n_after = adata.n_obs

    with open(output / "doublet_summary.tsv", "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"cells_before\t{n_before}\n")
        f.write(f"doublets_removed\t{n_before - n_after}\n")
        f.write(f"cells_after\t{n_after}\n")

    adata.write_h5ad(
        output / "doublet_filtered.h5ad"
    )


if __name__ == "__main__":
    main()