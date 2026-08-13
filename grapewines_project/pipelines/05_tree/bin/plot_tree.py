#!/usr/bin/env python3
"""
Draw the rooted, bootstrapped tree, coloured by the K = 7 assignment.

Two figures, because they answer different questions:

    fan          all 413 accessions at once, to see whether the K groups
                 form coherent clades. No tip labels: 413 names around a
                 circle hide the pattern the figure exists to show.

    rectangular  every tip labelled and every branch length honest, for
                 looking up individual accessions. Tall on purpose.

The tree is re-rooted on the outgroup here as well, even though SNPhylo
was given -o. It is a cheap safeguard, and phangorn's bootstrap step
returns whatever rooting its own tree object happened to carry.

Output:
    tree_<layout>.pdf
    tree_<layout>.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import newick


#
# The palette the task asks for, and the one the ADMIXTURE barplots
# already use: blue, orange, green, yellow, red, purple, cyan.
#

K_COLOURS = {
    "K1": "#00699A",
    "K2": "#FD702B",
    "K3": "#06402B",
    "K4": "#FDC700",
    "K5": "#AE0039",
    "K6": "#6A005C",
    "K7": "#08BDBD",
}

ADMIXED_COLOUR = "#9E9E9E"

OUTGROUP_COLOUR = "#000000"

BRANCH_COLOUR = "#666666"


def parse_args():

    parser = argparse.ArgumentParser(
        description="Draw the rooted bootstrap tree"
    )

    parser.add_argument("--tree", required=True, type=Path)

    parser.add_argument("--annotation", required=True, type=Path)

    parser.add_argument("--labels", type=Path, default=None)

    parser.add_argument(
        "--layout",
        choices=["fan", "rectangular"],
        default="fan",
    )

    parser.add_argument("--outgroup", default="ZZ01")

    parser.add_argument(
        "--min-support",
        type=float,
        default=70.0,
        help="Bootstrap values below this are not drawn",
    )

    parser.add_argument(
        "--drop-admixed",
        action="store_true",
        help="Prune tips below the Q threshold, leaving the seven "
             "groups and the outgroup only",
    )

    parser.add_argument("--outdir", required=True, type=Path)

    return parser.parse_args()


def normalise(name: str) -> str:

    first, sep, second = name.partition("_")

    return first if sep and first == second else name


def caption_subject(n_tips: int, args) -> str:
    """
    What the figure is actually showing -- which is not always all 413
    accessions, so the title has to say which.
    """

    if args.drop_admixed:

        return (
            f"Rooted ML phylogeny, {n_tips - 1} assigned accessions "
            f"+ {args.outgroup} (admixed samples pruned)"
        )

    return (
        f"Rooted ML phylogeny of {n_tips - 1} Caucasian accessions "
        f"+ {args.outgroup}"
    )


# --------------------------------------------------------------------- #
#  Layout
# --------------------------------------------------------------------- #

def assign_coordinates(root: newick.Node) -> int:
    """
    y spreads the tips evenly, x is the distance from the root.
    """

    sys.setrecursionlimit(200_000)

    counter = {"next": 0}

    def spread(node):

        if node.is_tip:

            node.y = float(counter["next"])

            counter["next"] += 1

            return node.y

        node.y = float(
            np.mean([spread(child) for child in node.children])
        )

        return node.y

    spread(root)

    def depth(node, so_far):

        node.x = so_far + node.length

        for child in node.children:
            depth(child, node.x)

    depth(root, 0.0)

    return counter["next"]


def draw_rectangular(root, n_tips, colours, args):

    height = max(12.0, 0.30 * n_tips)

    figure, axes = plt.subplots(figsize=(18, height))

    for node in newick.all_nodes(root):

        for child in node.children:

            axes.plot(
                [node.x, node.x], [node.y, child.y],
                color=BRANCH_COLOUR, linewidth=1.1, solid_capstyle="round",
            )

            axes.plot(
                [node.x, child.x], [child.y, child.y],
                color=BRANCH_COLOUR, linewidth=1.1, solid_capstyle="round",
            )

    tips = newick.tips(root)

    span = max(tip.x for tip in tips)

    axes.scatter(
        [tip.x for tip in tips],
        [tip.y for tip in tips],
        s=34,
        c=[colours[normalise(tip.name)] for tip in tips],
        zorder=3,
        linewidths=0,
    )

    for tip in tips:

        axes.text(
            tip.x + span * 0.004,
            tip.y,
            tip.name,
            fontsize=8.5,
            va="center",
            color=colours[normalise(tip.name)],
        )

    #
    # Bootstrap support, only where it is worth reading. Below 70 the
    # value should not be leaned on, and drawing all 411 of them turns
    # the figure into noise.
    #

    drawn = 0

    for node in newick.all_nodes(root):

        if node.is_tip or node.support is None:
            continue

        if node.support < args.min_support:
            continue

        axes.text(
            node.x - span * 0.003,
            node.y,
            f"{node.support:.0f}",
            fontsize=7.5,
            ha="right",
            va="bottom",
            color="#333333",
        )

        drawn += 1

    print(f"{drawn} bootstrap values at or above {args.min_support:.0f} drawn")

    axes.set_yticks([])

    axes.set_xlabel(
        "Substitutions per site (maximum likelihood)",
        fontsize=20,
    )

    axes.tick_params(axis="x", labelsize=16)

    axes.set_title(
        f"{caption_subject(n_tips, args)}\n"
        f"tips coloured by K = 7 ancestry; bootstrap shown where "
        f">= {args.min_support:.0f} of 100"
        + (
            "; support is from the full 413-tip tree"
            if args.drop_admixed else ""
        ),
        fontsize=24,
    )

    for side in ("top", "right", "left"):
        axes.spines[side].set_visible(False)

    return figure, axes


def compress_root_stem(root, outgroup: str, keep: float = 0.25):
    """
    ZZ01 is a different species, so the branch separating it from the
    Caucasian accessions is long -- here it is longer than the whole
    ingroup is deep.

    Rooting halves that branch and puts one half under each side, so the
    cost is not that ZZ01 sits far out: it is that every ingroup tip is
    pushed outwards by the same large constant. Drawn to scale the 412
    accessions end up in a thin ring at the rim with an empty disc in
    the middle, and none of the structure is legible.

    Both root branches are therefore cut back to `keep` times the depth
    of the ingroup itself, and the whole ingroup slides inwards with
    them. Relative distances *within* the ingroup, which is what the
    figure is read for, are untouched. Returns True when this happened
    so the caption can say so -- a shortened branch that is not labelled
    as shortened is a lie about the distances.
    """

    if len(root.children) != 2:
        return False

    outgroup_side = None
    ingroup_side = None

    for child in root.children:

        names = {normalise(name) for name in newick.tip_names(child)}

        if names == {outgroup}:
            outgroup_side = child

        else:
            ingroup_side = child

    if outgroup_side is None or ingroup_side is None:
        return False

    base = ingroup_side.x

    ingroup_nodes = newick.all_nodes(ingroup_side)

    spread = max(node.x for node in ingroup_nodes) - base

    if spread <= 0 or base <= spread * 0.5:
        return False

    delta = base - spread * keep

    for node in ingroup_nodes:
        node.x -= delta

    for node in newick.all_nodes(outgroup_side):
        node.x -= delta

    return True


def draw_fan(root, n_tips, colours, args):

    figure, axes = plt.subplots(figsize=(18, 18))

    # A wedge left open so the root direction is visible.
    span = 2 * np.pi * 0.96

    def angle(node):
        return node.y / max(n_tips - 1, 1) * span

    for node in newick.all_nodes(root):

        if node.is_tip:
            continue

        child_angles = [angle(child) for child in node.children]

        arc = np.linspace(min(child_angles), max(child_angles), 60)

        axes.plot(
            node.x * np.cos(arc),
            node.x * np.sin(arc),
            color=BRANCH_COLOUR,
            linewidth=0.9,
        )

        for child in node.children:

            theta = angle(child)

            axes.plot(
                [node.x * np.cos(theta), child.x * np.cos(theta)],
                [node.x * np.sin(theta), child.x * np.sin(theta)],
                color=BRANCH_COLOUR,
                linewidth=0.4,
            )

    tips = newick.tips(root)

    positions = np.array(
        [[tip.x * np.cos(angle(tip)), tip.x * np.sin(angle(tip))]
         for tip in tips]
    )

    axes.scatter(
        positions[:, 0],
        positions[:, 1],
        s=55,
        c=[colours[normalise(tip.name)] for tip in tips],
        zorder=3,
        linewidths=0,
    )

    #
    # The outgroup is the whole point of this figure being rooted, so it
    # gets a marker and a label of its own rather than one dot among 413.
    #

    radius = max(tip.x for tip in tips)

    for index, tip in enumerate(tips):

        if normalise(tip.name) != args.outgroup:
            continue

        x, y = positions[index]

        axes.scatter(
            [x], [y],
            s=380, marker="*", color=OUTGROUP_COLOUR, zorder=4, linewidths=0,
        )

        theta = angle(tip)

        axes.annotate(
            f"{tip.name}  (outgroup)",
            xy=(x, y),
            xytext=(
                (radius * 1.18) * np.cos(theta),
                (radius * 1.18) * np.sin(theta),
            ),
            fontsize=24,
            fontweight="bold",
            ha="center",
            va="center",
            arrowprops={
                "arrowstyle": "-",
                "color": OUTGROUP_COLOUR,
                "linewidth": 1.8,
            },
        )

    axes.set_aspect("equal")

    axes.axis("off")

    subtitle = "tips coloured by K = 7 ancestry (max Q >= 0.75)"

    if getattr(args, "outgroup_shortened", False):

        subtitle += "; the ZZ01 branch is drawn shortened"

    axes.set_title(
        f"{caption_subject(n_tips, args)}, rooted on {args.outgroup}\n"
        f"{subtitle}",
        fontsize=28,
    )

    return figure, axes


# --------------------------------------------------------------------- #

def main():

    args = parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)


    tree = newick.parse(args.tree.read_text())

    names = [normalise(name) for name in newick.tip_names(tree)]

    if args.outgroup not in names:
        raise ValueError(
            f"'{args.outgroup}' is not a tip of {args.tree.name}; refusing "
            f"to draw an unrooted tree as a rooted one"
        )

    tree = newick.reroot_on_tip(tree, args.outgroup)

    if newick.rescale_support(tree):
        print("support values were proportions; rescaled to percent")


    annotation = pd.read_csv(args.annotation, sep="\t", dtype={"tip": str})

    category = dict(zip(annotation["tip"], annotation["category"]))

    missing = [name for name in names if name not in category]

    if missing:
        raise ValueError(
            f"{len(missing)} tips have no annotation: "
            f"{', '.join(missing[:10])}"
        )


    #
    # Optionally reduce the figure to the seven assigned groups plus the
    # outgroup. Pruning after re-rooting, so the root stays on ZZ01.
    #
    # The support values are not recomputed: the bootstrap ran on all
    # 413 tips, so every number still describes a split of the full
    # tree. The caption says so.
    #

    if args.drop_admixed:

        drop = {
            name for name in names
            if category.get(name) == "Admixed"
        }

        tree = newick.prune_tips(tree, drop)

        names = [normalise(name) for name in newick.tip_names(tree)]

        print(f"dropped {len(drop)} admixed tips, {len(names)} left")


    n_tips = assign_coordinates(tree)

    print(f"{n_tips} tips, rooted on {args.outgroup}")


    def colour_of(name: str) -> str:

        value = category[name]

        if value == "Outgroup":
            return OUTGROUP_COLOUR

        if value == "Admixed":
            return ADMIXED_COLOUR

        return K_COLOURS.get(value, ADMIXED_COLOUR)

    colours = {name: colour_of(name) for name in names}


    #
    # Legend text. The K groups are described from the metadata, not
    # from a fixed table of K numbers: the column order of a .Q file is
    # arbitrary, so K3 means whatever this particular run made it mean.
    #

    descriptions = {}

    if args.labels is not None and args.labels.exists():

        labels = pd.read_csv(args.labels, sep="\t")

        descriptions = dict(zip(labels["group"], labels["label"]))

        counts = dict(zip(labels["group"], labels["n"]))

    else:

        counts = annotation["category"].value_counts().to_dict()


    def legend_text(group: str) -> str:

        n = counts.get(group, 0)

        description = descriptions.get(group, "")

        return f"{group} - {description}  (n = {n})" if description \
            else f"{group}  (n = {n})"

    #
    # Only categories that are actually on the figure.
    #

    drawn = {category[name] for name in names}

    entries = [
        (group, K_COLOURS[group])
        for group in sorted(K_COLOURS)
        if group in drawn
    ]

    if "Admixed" in drawn:
        entries.append(("Admixed", ADMIXED_COLOUR))

    if "Outgroup" in drawn:
        entries.append(("Outgroup", OUTGROUP_COLOUR))


    #
    # Only the fan gets the outgroup branch shortened. The rectangular
    # figure is the one that has to keep branch lengths honest, and it
    # has room to be as tall as the distances require.
    #

    args.outgroup_shortened = (
        compress_root_stem(tree, args.outgroup)
        if args.layout == "fan" else False
    )

    if args.outgroup_shortened:
        print("root stem compressed for the fan layout")

    draw = draw_fan if args.layout == "fan" else draw_rectangular

    figure, axes = draw(tree, n_tips, colours, args)

    axes.legend(
        handles=[
            Line2D(
                [], [],
                marker="*" if group == "Outgroup" else "o",
                linestyle="",
                markersize=20 if group == "Outgroup" else 16,
                color=colour,
                label=legend_text(group),
            )
            for group, colour in entries
        ],
        loc="upper left",
        # Outside the axes in both layouts: the fan fills its square to
        # the corners often enough that an inset legend lands on tips.
        bbox_to_anchor=(1.0, 1.0),
        frameon=False,
        fontsize=19,
        title="K = 7 ADMIXTURE assignment",
        title_fontsize=21,
    )

    figure.tight_layout()


    dpi = 300 if args.layout == "fan" else 150

    suffix = "_groups_only" if args.drop_admixed else ""

    pdf_path = args.outdir / f"tree_{args.layout}{suffix}.pdf"

    png_path = args.outdir / f"tree_{args.layout}{suffix}.png"

    figure.savefig(pdf_path, bbox_inches="tight")

    figure.savefig(png_path, dpi=dpi, bbox_inches="tight")

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
