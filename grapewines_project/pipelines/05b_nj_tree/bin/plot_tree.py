#!/usr/bin/env python3
"""
Draw a Newick tree, with tips coloured by metadata.

With 412 samples a rectangular phylogram would be about a metre tall, so
the default layout is circular. Rectangular is available for subsets.

Tip labels are off by default: 412 names around a circle are unreadable
and hide the structure the figure is meant to show. The colour of each
tip, and the legend, carry the information instead.

Input:
    tree.nwk
    metadata CSV, and/or sample_q_values_K<K>.tsv from 03_fst

Output:
    tree_<colour column>.pdf / .png
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


#
# The project palette, so a K group is the same colour in the tree as in
# the ADMIXTURE barplots and the FST figures.
#

CUSTOM_PALETTE = [
    "#00699A",
    "#FD702B",
    "#06402B",
    "#FDC700",
    "#AE0039",
    "#6A005C",
    "#08BDBD",
    "#F21B3F",
]

FALLBACK_PALETTE = [
    "#4E6E80", "#B07D62", "#7A9E7E", "#C4A24C",
    "#8C5F73", "#5B6C9E", "#A0522D", "#6B8E23",
    "#9DB4C0", "#D08770",
]

UNKNOWN_COLOUR = "#CFCFCF"


ID_COLUMN_CANDIDATES = [
    "Column 1", "IID", "ID", "id", "sample", "Sample", "tip"
]


# --------------------------------------------------------------------- #
#  Newick
# --------------------------------------------------------------------- #

class Node:

    __slots__ = ("name", "length", "support", "children", "x", "y")

    def __init__(self):
        self.name = ""
        self.length = 0.0
        self.support = None
        self.children = []
        self.x = 0.0
        self.y = 0.0


TOKEN = re.compile(r"[(),;]|[^(),;:]+|:[^(),;]+")


def parse_newick(text: str) -> Node:
    """
    Small recursive-descent parser. Handles branch lengths and the
    numeric internal labels that build_tree.py writes as support.
    """

    text = text.strip()

    position = 0

    def parse_node():

        nonlocal position

        node = Node()

        if text[position] == "(":

            position += 1

            while True:

                node.children.append(parse_node())

                if text[position] == ",":
                    position += 1
                    continue

                if text[position] == ")":
                    position += 1
                    break

                raise ValueError(
                    f"Unexpected '{text[position]}' at offset {position}"
                )

        # Label: a tip name, or a support value on an internal node.
        start = position

        while position < len(text) and text[position] not in "(),:;":
            position += 1

        label = text[start:position].strip()

        if node.children and label:

            try:
                node.support = float(label)
            except ValueError:
                node.name = label

        else:
            node.name = label

        if position < len(text) and text[position] == ":":

            position += 1

            start = position

            while position < len(text) and text[position] not in "(),;":
                position += 1

            node.length = float(text[start:position])

        return node

    root = parse_node()

    return root


def tips_of(node: Node) -> list[Node]:

    if not node.children:
        return [node]

    found = []

    for child in node.children:
        found.extend(tips_of(child))

    return found


def all_nodes(node: Node) -> list[Node]:

    found = [node]

    for child in node.children:
        found.extend(all_nodes(child))

    return found


# --------------------------------------------------------------------- #
#  Layout
# --------------------------------------------------------------------- #

def assign_coordinates(root: Node, scale: str = "linear"):
    """
    x is the horizontal position, y spreads the tips evenly and puts
    each internal node at the mean of its children.

    scale="linear"
        x is the true distance from the root. Faithful, but on an
        individual-level NJ tree the private branch of each sample is
        far longer than the splits between groups, so everything
        interesting collapses towards the centre.

    scale="cladogram"
        x is the node's level, with every tip pushed to the same depth.
        Branch lengths are discarded and only the topology remains.
    """

    sys.setrecursionlimit(100_000)

    counter = {"next": 0}

    def spread(node):

        if not node.children:

            node.y = counter["next"]

            counter["next"] += 1

            return node.y

        positions = [spread(child) for child in node.children]

        node.y = float(np.mean(positions))

        return node.y

    spread(root)


    if scale == "cladogram":

        def height_of(node):

            if not node.children:
                node.x = 0.0
                return 0

            height = 1 + max(height_of(child) for child in node.children)

            node.x = float(height)

            return height

        total = height_of(root)

        for node in all_nodes(root):
            node.x = total - node.x

    else:

        def depth_of(node, depth):

            node.x = depth + node.length

            for child in node.children:
                depth_of(child, node.x)

        depth_of(root, 0.0)

    return counter["next"]


def draw_rectangular(root, n_tips, colours, labels, show_labels, title, x_label=True):

    height = max(4.0, 0.14 * n_tips)

    figure, axes = plt.subplots(figsize=(10, height))

    for node in all_nodes(root):

        for child in node.children:

            # Elbow: vertical at the parent's x, then horizontal out.
            axes.plot(
                [node.x, node.x],
                [node.y, child.y],
                color="#666666",
                linewidth=0.6
            )

            axes.plot(
                [node.x, child.x],
                [child.y, child.y],
                color="#666666",
                linewidth=0.6
            )

    tips = tips_of(root)

    axes.scatter(
        [tip.x for tip in tips],
        [tip.y for tip in tips],
        s=12,
        c=[colours[tip.name] for tip in tips],
        zorder=3,
        linewidths=0
    )

    if show_labels:

        span = max(tip.x for tip in tips)

        for tip in tips:

            axes.text(
                tip.x + span * 0.005,
                tip.y,
                tip.name,
                fontsize=5,
                va="center"
            )

    axes.set_yticks([])

    axes.set_xlabel(
        "Genetic distance (1 - IBS)"
        if x_label else "Topology only, branch lengths discarded"
    )

    axes.set_title(title, fontsize=11)

    for side in ("top", "right", "left"):
        axes.spines[side].set_visible(False)

    return figure, axes


def draw_circular(root, n_tips, colours, labels, show_labels, title, x_label=True):

    figure, axes = plt.subplots(figsize=(11, 11))

    # A wedge is left open so the tree reads as a fan rather than a
    # closed ring, which makes the root direction obvious.
    span = 2 * np.pi * 0.97

    def angle_of(node):
        return node.y / max(n_tips - 1, 1) * span

    def to_xy(node):

        theta = angle_of(node)

        return node.x * np.cos(theta), node.x * np.sin(theta)


    for node in all_nodes(root):

        if not node.children:
            continue

        # Arc joining the children, drawn at the parent's radius.
        child_angles = [angle_of(child) for child in node.children]

        arc = np.linspace(min(child_angles), max(child_angles), 40)

        axes.plot(
            node.x * np.cos(arc),
            node.x * np.sin(arc),
            color="#666666",
            linewidth=0.5
        )

        for child in node.children:

            theta = angle_of(child)

            axes.plot(
                [node.x * np.cos(theta), child.x * np.cos(theta)],
                [node.x * np.sin(theta), child.x * np.sin(theta)],
                color="#666666",
                linewidth=0.5
            )


    tips = tips_of(root)

    positions = np.array([to_xy(tip) for tip in tips])

    axes.scatter(
        positions[:, 0],
        positions[:, 1],
        s=14,
        c=[colours[tip.name] for tip in tips],
        zorder=3,
        linewidths=0
    )

    if show_labels:

        radius = max(tip.x for tip in tips)

        for tip in tips:

            theta = angle_of(tip)

            degrees = np.degrees(theta)

            axes.text(
                (tip.x + radius * 0.02) * np.cos(theta),
                (tip.x + radius * 0.02) * np.sin(theta),
                tip.name,
                fontsize=4,
                rotation=degrees if -90 < degrees < 90 else degrees + 180,
                rotation_mode="anchor",
                ha="left" if -90 < degrees < 90 else "right",
                va="center"
            )

    axes.set_aspect("equal")

    axes.axis("off")

    axes.set_title(title, fontsize=12)

    return figure, axes


# --------------------------------------------------------------------- #

def parse_args():

    parser = argparse.ArgumentParser(
        description="Draw a Newick tree coloured by metadata"
    )

    parser.add_argument("--tree", required=True, type=Path)

    parser.add_argument(
        "--metadata",
        type=Path,
        default=None,
        help="Sample metadata CSV"
    )

    parser.add_argument(
        "--groups",
        type=Path,
        default=None,
        help="sample_q_values_K<K>.tsv from 03_fst, to colour by "
             "ADMIXTURE group"
    )

    parser.add_argument(
        "--color-by",
        default="Geographic origin by country",
        help="Metadata column, or 'assignment' for the ADMIXTURE group"
    )

    parser.add_argument(
        "--layout",
        choices=["circular", "rectangular"],
        default="circular"
    )

    parser.add_argument(
        "--scale",
        choices=["linear", "cladogram"],
        default="linear",
        help="linear keeps branch lengths; cladogram drops them and "
             "shows the topology only"
    )

    parser.add_argument(
        "--show-labels",
        action="store_true",
        help="Draw tip names; unreadable above about 80 samples"
    )

    parser.add_argument("--outdir", required=True, type=Path)

    return parser.parse_args()


def find_id_column(table: pd.DataFrame) -> str:

    for column in ID_COLUMN_CANDIDATES:
        if column in table.columns:
            return column

    raise ValueError(
        f"No sample ID column; available: {list(table.columns)}"
    )


def build_colour_map(tips, annotation, column):
    """
    Returns (tip -> colour, level -> colour) with a stable level order.
    """

    values = {
        tip: str(annotation.get(tip, "")).strip()
        for tip in tips
    }

    levels = sorted(
        {
            value for value in values.values()
            if value and value.upper() not in ("NA", "NAN")
        }
    )

    #
    # ADMIXTURE groups keep the project palette in their own order, so
    # K3 in the tree is the K3 colour everywhere else.
    #

    if all(re.fullmatch(r"K\d+", level) for level in levels) and levels:

        levels = sorted(levels, key=lambda s: int(s[1:]))

        palette = CUSTOM_PALETTE

    else:

        palette = (
            CUSTOM_PALETTE
            if len(levels) <= len(CUSTOM_PALETTE)
            else FALLBACK_PALETTE
        )

    legend = {
        level: palette[i % len(palette)]
        for i, level in enumerate(levels)
    }

    colours = {
        tip: legend.get(value, UNKNOWN_COLOUR)
        for tip, value in values.items()
    }

    if any(colour == UNKNOWN_COLOUR for colour in colours.values()):
        legend["no data"] = UNKNOWN_COLOUR

    return colours, legend


def main():

    args = parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)


    root = parse_newick(args.tree.read_text())

    n_tips = assign_coordinates(root, args.scale)

    tips = [tip.name for tip in tips_of(root)]

    print(f"{n_tips} tips")


    #
    # Where the colouring comes from.
    #

    annotation = {}

    if args.color_by == "assignment":

        if args.groups is None:
            raise ValueError(
                "--color-by assignment needs --groups"
            )

        table = pd.read_csv(args.groups, sep="\t", dtype={"IID": str})

        annotation = dict(zip(table["IID"], table["assignment"]))

    elif args.metadata is not None:

        table = pd.read_csv(args.metadata)

        id_column = find_id_column(table)

        if args.color_by not in table.columns:
            raise ValueError(
                f"Column '{args.color_by}' is not in the metadata. "
                f"Available: {list(table.columns)}"
            )

        annotation = dict(
            zip(
                table[id_column].astype(str),
                table[args.color_by]
            )
        )

    matched = sum(1 for tip in tips if tip in annotation)

    print(f"{matched} of {len(tips)} tips matched to '{args.color_by}'")

    if matched == 0:
        raise ValueError(
            "No tip name matches the annotation table; check that the "
            ".fam IDs and the metadata IDs are the same"
        )


    colours, legend = build_colour_map(tips, annotation, args.color_by)


    title = (
        f"Neighbour-joining tree, coloured by {args.color_by}"
        + ("  (cladogram: branch lengths not to scale)"
           if args.scale == "cladogram" else "")
    )

    draw = (
        draw_circular
        if args.layout == "circular"
        else draw_rectangular
    )

    figure, axes = draw(
        root,
        n_tips,
        colours,
        legend,
        args.show_labels,
        title,
        args.scale == "linear"
    )

    axes.legend(
        handles=[
            Line2D(
                [], [],
                marker="o",
                linestyle="",
                markersize=6,
                color=colour,
                label=level
            )
            for level, colour in legend.items()
        ],
        loc="upper left",
        bbox_to_anchor=(1.0, 1.0) if args.layout == "rectangular" else (0.98, 1.0),
        frameon=False,
        fontsize=8
    )

    figure.tight_layout()


    safe = re.sub(r"[^A-Za-z0-9]+", "_", args.color_by).strip("_")

    if args.scale == "cladogram":
        safe += "_cladogram"

    pdf_path = args.outdir / f"tree_{safe}.pdf"
    png_path = args.outdir / f"tree_{safe}.png"

    figure.savefig(pdf_path, bbox_inches="tight")
    figure.savefig(png_path, dpi=200, bbox_inches="tight")

    plt.close(figure)

    print("Wrote:")
    print(" ", pdf_path)
    print(" ", png_path)


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
