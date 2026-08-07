#!/usr/bin/env python3

import argparse
from pathlib import Path

import scanpy as sc
import harmonypy


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--n-pcs", type=int, default=30)

    args = parser.parse_args()

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.input)

    if "condition" not in adata.obs:
        raise ValueError(
            "Column 'condition' is required for Harmony"
        )

    if "X_pca" not in adata.obsm:
        raise ValueError(
            "PCA is missing. Run normalize.py first."
        )

    pca = adata.obsm["X_pca"][:, :args.n_pcs]

    harmony_result = harmonypy.run_harmony(
        pca,
        adata.obs,
        vars_use=["condition"],
    )

    Z_corr = harmony_result.Z_corr

    # harmonypy versions differ:
    # some return cells x PCs, others PCs x cells

    if Z_corr.shape[0] == adata.n_obs:
        adata.obsm["X_harmony"] = Z_corr

    elif Z_corr.shape[1] == adata.n_obs:
        adata.obsm["X_harmony"] = Z_corr.T

    else:
        raise ValueError(
            f"Unexpected Harmony output shape: {Z_corr.shape}. "
            f"Expected one dimension to match number of cells: {adata.n_obs}"
        )

    adata.write_h5ad(
        output / "harmony.h5ad"
    )


if __name__ == "__main__":
    main()