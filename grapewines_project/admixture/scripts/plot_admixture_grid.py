#!/usr/bin/env python3
"""
Plot ADMIXTURE results for every K, with a shared sample ordering
defined by an order file, and metadata colour tracks below.

Improvements
------------
1. Custom 8-colour palette matching the user's R palette.
2. Greedy component alignment across K values so that the "same"
   ancestry cluster keeps the same colour at every K.
3. Optional --cv flag to produce a cross-validation error plot
   (CV error vs K) to justify the chosen K.
4. fill_between(step='post') for fast, correct rendering.
5. rasterized=True for lightweight vector output.
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


# ------------------------------------------------------------------ #
#  Custom palette (matches R: my_palette)
# ------------------------------------------------------------------ #

CUSTOM_PALETTE = [
    "#00699A",  # K1 — teal-blue
    "#FD702B",  # K2 — orange
    "#06402B",  # K3 — dark green
    "#FDC700",  # K4 — gold
    "#AE0039",  # K5 — crimson
    "#6A005C",  # K6 — purple
    "#08BDBD",  # K7 — cyan
    "#F21B3F",  # K8 — red
]


# ------------------------------------------------------------------ #
#  CLI
# ------------------------------------------------------------------ #

def parse_args():

    p = argparse.ArgumentParser(
        description="Grid of ADMIXTURE barplots across K values"
    )

    p.add_argument("--qdir",     type=Path, required=True,
                   help="Directory with *.Q files")
    p.add_argument("--orders",   type=Path, required=True,
                   help="Directory with order_Q*.tsv files")
    p.add_argument("--fam",      type=Path, required=True,
                   help="PLINK .fam matching Q rows")
    p.add_argument("--metadata", type=Path, required=True,
                   help="CSV with sample metadata")
    p.add_argument("--outdir",   type=Path, required=True)
    p.add_argument("--cv",       type=Path, default=None,
                   help="Path to ADMIXTURE log file or CV-error TSV "
                        "(columns: K, CV) to produce CV-error plot")
    p.add_argument("--best-k",   type=int,  default=None,
                   help="Highlight this K on the CV-error plot "
                        "(default: K with minimum CV error)")

    return p.parse_args()


# ------------------------------------------------------------------ #
#  I/O helpers
# ------------------------------------------------------------------ #

def get_k(path):
    return int(re.search(r"\.(\d+)\.Q", path.name).group(1))


def read_q(path):
    return pd.read_csv(path, sep=r"\s+", header=None).values


def read_fam(path):
    fam = pd.read_csv(path, sep=r"\s+", header=None)
    fam.columns = ["FID", "IID", "PID", "MID", "SEX", "PHENO"]
    return fam


def read_order(path):
    return pd.read_csv(path, sep="\t")["IID"].astype(str).tolist()


def find_metadata_id_column(metadata):
    for col in ("Column 1", "IID", "ID", "id", "sample", "Sample"):
        if col in metadata.columns:
            return col
    raise ValueError(
        f"Cannot find sample ID column. "
        f"Available: {list(metadata.columns)}"
    )


# ------------------------------------------------------------------ #
#  CV-error parsing
# ------------------------------------------------------------------ #

def parse_cv_errors(path):
    """
    Accept either:
      • An ADMIXTURE stdout log containing lines like
            CV error (K=3): 0.42857
      • A two-column TSV / CSV with header K,CV
    Returns a DataFrame with columns [K, CV] sorted by K.
    """
    text = path.read_text()

    # Try ADMIXTURE log format first
    pattern = r"CV error \(K=(\d+)\):\s+([0-9.]+)"
    matches = re.findall(pattern, text)

    if matches:
        df = pd.DataFrame(matches, columns=["K", "CV"])
        df["K"]  = df["K"].astype(int)
        df["CV"] = df["CV"].astype(float)
    else:
        # Fall back to tabular format
        sep = "\t" if "\t" in text else ","
        df = pd.read_csv(path, sep=sep)
        # Normalise column names
        df.columns = [c.strip() for c in df.columns]
        if "K" not in df.columns or "CV" not in df.columns:
            raise ValueError(
                f"CV file must have columns K and CV. "
                f"Found: {list(df.columns)}"
            )
        df["K"]  = df["K"].astype(int)
        df["CV"] = df["CV"].astype(float)

    return df.sort_values("K").reset_index(drop=True)


# ------------------------------------------------------------------ #
#  Component alignment across K values
# ------------------------------------------------------------------ #

def align_components(q_frames, sorted_ks):
    """
    Greedy alignment of ancestry components across increasing K.

    For K and K-1, compute the N×(K-1) overlap matrix between every
    pair of components (sum of element-wise min across all samples).
    Greedily assign each new component to the colour of the K-1
    component it overlaps most with; the leftover component gets the
    next unused colour slot.

    Returns
    -------
    colour_orders : dict[int, list[int]]
        For each K, a list of length K giving the colour index
        (0-based into the palette) for each column of the Q matrix.
    """

    colour_orders = {}

    # K=2 is the base: assign colours 0, 1
    k0 = sorted_ks[0]
    colour_orders[k0] = list(range(k0))

    next_colour = k0  # next unused palette slot

    for idx in range(1, len(sorted_ks)):
        k_prev = sorted_ks[idx - 1]
        k_cur  = sorted_ks[idx]

        q_prev = q_frames[k_prev]  # (N, k_prev)
        q_cur  = q_frames[k_cur]   # (N, k_cur)

        prev_order = colour_orders[k_prev]

        # Overlap matrix: (k_cur, k_prev)
        # overlap[i, j] = sum_n min(q_cur[n, i], q_prev[n, j])
        overlap = np.zeros((k_cur, k_prev))
        for i in range(k_cur):
            for j in range(k_prev):
                overlap[i, j] = np.minimum(
                    q_cur[:, i], q_prev[:, j]
                ).sum()

        # Greedy matching: assign each cur component to the best
        # matching prev component (by overlap), breaking ties by index
        assigned_cur   = set()          # cur components already matched
        assigned_color = set()          # palette colours already used

        cur_colours = [None] * k_cur

        # Sort assignments by descending overlap
        pairs = []
        for i in range(k_cur):
            for j in range(k_prev):
                pairs.append((overlap[i, j], i, j))
        pairs.sort(key=lambda x: -x[0])

        for _, i, j in pairs:
            if i in assigned_cur:
                continue
            c = prev_order[j]
            if c in assigned_color:
                continue
            cur_colours[i] = c
            assigned_cur.add(i)
            assigned_color.add(c)

        # Any unassigned components get the next fresh colour
        for i in range(k_cur):
            if cur_colours[i] is None:
                cur_colours[i] = next_colour
                next_colour += 1

        colour_orders[k_cur] = cur_colours

    return colour_orders


# ------------------------------------------------------------------ #
#  Drawing primitives
# ------------------------------------------------------------------ #

def draw_admixture(ax, q, colour_indices, palette):
    """
    Stacked ancestry barplot via fill_between(step='post').

    colour_indices : list[int]
        Maps each Q column to a palette index so that the same
        ancestry component keeps the same colour across panels.
    """

    n, k = q.shape

    x = np.arange(n + 1)
    bottom = np.zeros(n + 1)

    for j in range(k):
        col = np.append(q[:, j], 0.0)
        top = bottom + col

        ax.fill_between(
            x, bottom, top,
            step="post",
            facecolor=palette[colour_indices[j] % len(palette)],
            edgecolor="none",
            linewidth=0,
            rasterized=True,
        )

        bottom = top

    ax.set_xlim(0, n)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xticks([])


def draw_metadata_track(ax, values, n_samples, palette):
    """Single horizontal colour bar for a categorical metadata column."""

    categories = sorted(values.unique().tolist())
    mapping = {v: i for i, v in enumerate(categories)}
    numeric = np.array([mapping[v] for v in values])

    cmap = ListedColormap(palette[: len(categories)])

    ax.imshow(
        [numeric],
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        extent=[0, n_samples, 0, 1],
        rasterized=True,
    )

    ax.set_yticks([])
    ax.set_xticks([])


# ------------------------------------------------------------------ #
#  CV-error plot
# ------------------------------------------------------------------ #

def plot_cv(cv_df, outdir, best_k=None):
    """
    Bar + line plot of cross-validation error vs K.

    Highlights the best K (minimum CV or user-specified) so that it is
    immediately obvious why that K was chosen.
    """

    if best_k is None:
        best_k = int(cv_df.loc[cv_df["CV"].idxmin(), "K"])

    ks   = cv_df["K"].values
    vals = cv_df["CV"].values

    fig, ax = plt.subplots(figsize=(7, 4.5))

    bar_colors = [
        CUSTOM_PALETTE[0] if k == best_k else "#BFBFBF"
        for k in ks
    ]

    bars = ax.bar(ks, vals, color=bar_colors, edgecolor="white",
                  linewidth=0.6, zorder=2)

    # Line overlay
    ax.plot(ks, vals, color="#333333", marker="o", markersize=5,
            linewidth=1.2, zorder=3)

    # Highlight best K
    best_val = float(cv_df.loc[cv_df["K"] == best_k, "CV"].iloc[0])
    ax.annotate(
        f"K = {best_k}\nCV = {best_val:.5f}",
        xy=(best_k, best_val),
        xytext=(best_k + 0.6, best_val + (vals.max() - vals.min()) * 0.12),
        fontsize=9,
        fontweight="bold",
        color=CUSTOM_PALETTE[0],
        arrowprops=dict(
            arrowstyle="->",
            color=CUSTOM_PALETTE[0],
            lw=1.3,
        ),
        zorder=4,
    )

    ax.set_xlabel("K  (number of ancestral populations)", fontsize=11)
    ax.set_ylabel("Cross-validation error", fontsize=11)
    ax.set_title("ADMIXTURE cross-validation error", fontsize=13,
                 fontweight="bold")

    ax.set_xticks(ks)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linewidth=0.3, alpha=0.5)

    fig.tight_layout()

    out = outdir / "cv_error.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"  saved {out}")


# ------------------------------------------------------------------ #
#  Main
# ------------------------------------------------------------------ #

def main():

    args = parse_args()
    args.outdir.mkdir(exist_ok=True, parents=True)

    # ---- FAM ----------------------------------------------------- #

    fam = read_fam(args.fam)
    fam_iids = fam["IID"].astype(str).tolist()

    # ---- Q matrices (indexed by FAM IIDs) ------------------------ #

    q_files = sorted(
        args.qdir.glob("cauc_filtered.final.*.Q"),
        key=get_k,
    )

    if not q_files:
        raise FileNotFoundError(f"No Q files found in {args.qdir}")

    q_frames = {}

    for f in q_files:
        k = get_k(f)
        arr = read_q(f)

        if len(arr) != len(fam_iids):
            raise ValueError(
                f"{f.name}: {len(arr)} rows vs "
                f"{len(fam_iids)} in FAM"
            )

        q_frames[k] = pd.DataFrame(arr, index=fam_iids)

    # ---- Align components across K values ------------------------ #

    sorted_ks = sorted(q_frames.keys())

    # We need raw numpy arrays for alignment
    q_arrays = {k: q_frames[k].values for k in sorted_ks}
    colour_orders = align_components(q_arrays, sorted_ks)

    # ---- Metadata ------------------------------------------------ #

    metadata = pd.read_csv(args.metadata)
    meta_id = find_metadata_id_column(metadata)
    metadata[meta_id] = metadata[meta_id].astype(str)
    metadata = metadata.set_index(meta_id)

    meta_cols = metadata.columns.tolist()
    n_meta = len(meta_cols)

    # ---- Palette ------------------------------------------------- #

    # Extend if more than 8 components are ever needed
    palette = list(CUSTOM_PALETTE)
    if max(max(v) for v in colour_orders.values()) >= len(palette):
        extra = plt.cm.get_cmap("Set3", 12)
        palette += [extra(i) for i in range(12)]

    meta_palette_cmap = plt.cm.get_cmap("tab20", 20)
    meta_palette = [meta_palette_cmap(i) for i in range(20)]

    # ---- CV-error plot ------------------------------------------- #

    if args.cv is not None:
        cv_df = parse_cv_errors(args.cv)
        plot_cv(cv_df, args.outdir, best_k=args.best_k)

    # ---- One figure per order file ------------------------------- #

    order_files = sorted(args.orders.glob("order_Q*.tsv"))

    if not order_files:
        raise FileNotFoundError(
            f"No order files found in {args.orders}"
        )

    for order_file in order_files:

        order_name = order_file.stem
        order = read_order(order_file)

        n_samples = len(order)
        n_q = len(q_files)
        n_rows = n_q + n_meta

        # ---- layout ---------------------------------------------- #

        q_h = 1.8
        meta_h = 0.35
        height_ratios = [q_h] * n_q + [meta_h] * n_meta

        width = max(14, n_samples / 12)
        height = sum(height_ratios) + 1.5

        fig, axes = plt.subplots(
            n_rows, 1,
            figsize=(width, height),
            gridspec_kw=dict(
                height_ratios=height_ratios,
                hspace=0.06,
            ),
        )

        if n_rows == 1:
            axes = [axes]

        # ---- ADMIXTURE rows -------------------------------------- #

        for i, qf in enumerate(q_files):
            k = get_k(qf)
            q_ordered = q_frames[k].loc[order].values

            draw_admixture(
                axes[i],
                q_ordered,
                colour_orders[k],
                palette,
            )

            axes[i].set_ylabel(
                f"K={k}",
                rotation=0, ha="right", va="center",
                fontsize=9,
            )
            axes[i].yaxis.set_label_coords(-0.02, 0.5)

        # ---- metadata rows --------------------------------------- #

        meta_ordered = metadata.reindex(order)

        for j, col in enumerate(meta_cols):
            ax_idx = n_q + j
            values = meta_ordered[col].fillna("NA").astype(str)
            draw_metadata_track(
                axes[ax_idx], values, n_samples, meta_palette
            )

            axes[ax_idx].set_ylabel(
                col,
                rotation=0, ha="right", va="center",
                fontsize=8,
            )
            axes[ax_idx].yaxis.set_label_coords(-0.02, 0.5)

        # ---- title & save ---------------------------------------- #

        ref_k = order_name.replace("order_", "")
        fig.suptitle(
            f"ADMIXTURE — sample ordering by {ref_k}",
            fontsize=13, y=0.995,
        )

        plt.subplots_adjust(
            left=0.10, right=0.98,
            top=0.97,  bottom=0.02,
        )

        out_png = args.outdir / f"{order_name}.png"
        fig.savefig(out_png, dpi=200, bbox_inches="tight")
        plt.close(fig)

        print(f"  saved {out_png}")

    print("Done.")


if __name__ == "__main__":
    main()