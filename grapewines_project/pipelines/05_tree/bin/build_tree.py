#!/usr/bin/env python3
"""
Neighbour-joining tree from SNP genotypes.

    .traw  ->  1-IBS distance matrix  ->  NJ  ->  midpoint root  ->  Newick

Distances are computed here rather than taken from PLINK so that the
bootstrap replicates use exactly the same estimator as the tree they
support. When --mdist is given, the matrix is checked against PLINK's
own `--distance 1-ibs` output, which catches a mistake in this code
rather than letting it through silently.

Bootstrap resamples loci with replacement, rebuilds the tree, and counts
how often each split of the reference tree reappears. Resampling is done
with multinomial weights instead of by copying columns, which keeps the
genotype matrix in memory once.

Output:
    tree.nwk                Newick, with bootstrap support on internal nodes
    distance_matrix.tsv     the 1-IBS matrix
    tree_tips.tsv           tip names, in the order used internally
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():

    parser = argparse.ArgumentParser(
        description="Neighbour-joining tree from a PLINK .traw"
    )

    parser.add_argument(
        "--traw",
        required=True,
        type=Path,
        help="PLINK --recode A-transpose output (SNPs in rows)"
    )

    parser.add_argument(
        "--mdist",
        type=Path,
        default=None,
        help="PLINK --distance square 1-ibs matrix, used as a cross-check"
    )

    parser.add_argument(
        "--bootstrap",
        type=int,
        default=100,
        help="Bootstrap replicates; 0 disables (default 100)"
    )

    parser.add_argument(
        "--root",
        choices=["midpoint", "outgroup", "none"],
        default="midpoint",
        help="How to root the tree (default midpoint)"
    )

    parser.add_argument(
        "--outgroup",
        default=None,
        help="Sample ID to root on, when --root outgroup"
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Seed for the bootstrap resampling"
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path
    )

    return parser.parse_args()


# --------------------------------------------------------------------- #
#  Genotypes and distances
# --------------------------------------------------------------------- #

def read_traw(path: Path) -> tuple[np.ndarray, list[str]]:
    """
    Returns (genotypes, sample names) with genotypes shaped
    (n_samples, n_loci) and missing calls as NaN.
    """

    table = pd.read_csv(path, sep="\t")

    #
    # The first six columns are CHR, SNP, (C)M, POS, COUNTED, ALT; the
    # rest are one column per sample, named FID_IID.
    #

    metadata_columns = 6

    if table.shape[1] <= metadata_columns:
        raise ValueError(
            f"{path} has {table.shape[1]} columns; expected the six "
            f"PLINK columns followed by one per sample"
        )

    sample_columns = list(table.columns[metadata_columns:])

    genotypes = (
        table[sample_columns]
        .to_numpy(dtype=np.float32)
        .T
        .copy()
    )

    #
    # PLINK writes FID_IID. FID equals IID everywhere in this dataset,
    # so the doubled name is dropped for readability.
    #

    names = []

    for column in sample_columns:

        fid, _, iid = column.partition("_")

        names.append(iid if iid and iid == fid else column)

    return genotypes, names


def one_hot(genotypes: np.ndarray):
    """
    Split the allele counts into indicator matrices, once, so every
    bootstrap replicate is four matrix products and nothing else.
    """

    present = ~np.isnan(genotypes)

    g0 = ((genotypes == 0) & present).astype(np.float32)
    g1 = ((genotypes == 1) & present).astype(np.float32)
    g2 = ((genotypes == 2) & present).astype(np.float32)

    return g0, g1, g2, present.astype(np.float32)


def ibs_distance(g0, g1, g2, present, weights=None) -> np.ndarray:
    """
    1 - IBS, averaged over the loci where both samples are called.

    For allele counts in {0, 1, 2} the per-locus contribution is
    |g_i - g_j| / 2, so the sum over loci decomposes into counts of the
    (0,1), (1,2) and (0,2) pairings. Each of those is one matrix
    product, and its transpose supplies the mirrored pairing.
    """

    if weights is None:

        w0, w1, w2, wm = g0, g1, g2, present

    else:

        w0 = g0 * weights
        w1 = g1 * weights
        w2 = g2 * weights
        wm = present * weights

    a01 = w0 @ g1.T
    a12 = w1 @ g2.T
    a02 = w0 @ g2.T

    difference = (
        (a01 + a01.T)
        + (a12 + a12.T)
        + 2.0 * (a02 + a02.T)
    )

    comparable = wm @ present.T

    with np.errstate(invalid="ignore", divide="ignore"):

        distance = difference / (2.0 * comparable)

    distance[comparable == 0] = np.nan

    np.fill_diagonal(distance, 0.0)

    return distance


# --------------------------------------------------------------------- #
#  Neighbour joining
# --------------------------------------------------------------------- #

def neighbour_joining(distance: np.ndarray) -> tuple[list, int]:
    """
    Saitou and Nei's algorithm.

    Returns (edges, n_nodes), where edges are (u, v, length) on an
    unrooted tree whose first n_tips nodes are the tips.
    """

    n_tips = distance.shape[0]

    if n_tips < 3:
        raise ValueError("Need at least three samples to build a tree")

    #
    # Joining n tips creates n - 2 internal nodes, so the working
    # matrix is sized for all of them up front.
    #

    n_nodes = 2 * n_tips - 2

    d = np.zeros((n_nodes, n_nodes))

    d[:n_tips, :n_tips] = distance

    active = list(range(n_tips))

    edges = []

    next_node = n_tips


    while len(active) > 2:

        index = np.array(active)

        sub = d[np.ix_(index, index)]

        m = len(active)

        totals = sub.sum(axis=1)

        #
        # Q_ij = (m - 2) d_ij - R_i - R_j, minimised over i != j.
        #

        q = (m - 2) * sub - totals[:, None] - totals[None, :]

        np.fill_diagonal(q, np.inf)

        flat = np.argmin(q)

        i, j = np.unravel_index(flat, q.shape)

        node_i = active[i]
        node_j = active[j]

        d_ij = sub[i, j]

        delta = (totals[i] - totals[j]) / (m - 2)

        limb_i = 0.5 * d_ij + 0.5 * delta
        limb_j = d_ij - limb_i

        # Negative limbs are an artefact of non-additive distances.
        limb_i = max(limb_i, 0.0)
        limb_j = max(limb_j, 0.0)

        edges.append((next_node, node_i, limb_i))
        edges.append((next_node, node_j, limb_j))

        #
        # Distance from the new node to every remaining cluster.
        #

        keep = np.ones(m, dtype=bool)

        keep[i] = False
        keep[j] = False

        remaining = index[keep]

        new_distances = 0.5 * (
            sub[i, keep] + sub[j, keep] - d_ij
        )

        d[next_node, remaining] = new_distances
        d[remaining, next_node] = new_distances

        d[next_node, next_node] = 0.0

        active = list(remaining) + [next_node]

        next_node += 1


    last_a, last_b = active

    edges.append((last_a, last_b, max(d[last_a, last_b], 0.0)))

    return edges, next_node


def adjacency(edges, n_nodes) -> list[list[tuple[int, float]]]:

    neighbours = [[] for _ in range(n_nodes)]

    for u, v, length in edges:

        neighbours[u].append((v, length))
        neighbours[v].append((u, length))

    return neighbours


# --------------------------------------------------------------------- #
#  Splits, for bootstrap support
# --------------------------------------------------------------------- #

def splits_of(edges, n_nodes, n_tips) -> dict[frozenset, tuple]:
    """
    Every internal edge of an unrooted tree cuts the tips into two sets.
    The split is named by whichever side leaves out tip 0, so the same
    bipartition gets the same name in any tree over the same tips.
    """

    neighbours = adjacency(edges, n_nodes)

    result = {}

    for u, v, length in edges:

        if u < n_tips or v < n_tips:
            # A pendant edge; its split is a single tip and carries no
            # topological information.
            continue

        # Tips reachable from v without crossing the u-v edge.
        seen = set()

        stack = [(v, u)]

        while stack:

            node, came_from = stack.pop()

            if node < n_tips:
                seen.add(node)

            for neighbour, _ in neighbours[node]:

                if neighbour != came_from:
                    stack.append((neighbour, node))

        side = frozenset(seen)

        if 0 in side:
            side = frozenset(set(range(n_tips)) - side)

        result[side] = (u, v, length)

    return result


# --------------------------------------------------------------------- #
#  Rooting
# --------------------------------------------------------------------- #

def path_distances_from(start, neighbours, n_nodes):

    distance = np.full(n_nodes, np.inf)

    parent = np.full(n_nodes, -1, dtype=int)

    distance[start] = 0.0

    stack = [start]

    while stack:

        node = stack.pop()

        for neighbour, length in neighbours[node]:

            if distance[neighbour] == np.inf:

                distance[neighbour] = distance[node] + length
                parent[neighbour] = node

                stack.append(neighbour)

    return distance, parent


def midpoint_of_tree(edges, n_nodes, n_tips):
    """
    The two most distant tips, and the point halfway along the path
    between them. Returns (edge_u, edge_v, distance_from_u).
    """

    neighbours = adjacency(edges, n_nodes)

    # Farthest tip from an arbitrary tip, then farthest from that one:
    # the usual two-pass diameter search, exact on a tree.

    first, _ = path_distances_from(0, neighbours, n_nodes)

    end_a = int(np.argmax(np.where(np.arange(n_nodes) < n_tips, first, -1)))

    from_a, parent = path_distances_from(end_a, neighbours, n_nodes)

    end_b = int(np.argmax(np.where(np.arange(n_nodes) < n_tips, from_a, -1)))

    diameter = from_a[end_b]

    half = diameter / 2.0

    # Walk back from end_b towards end_a until half the diameter is left
    # behind.

    node = end_b

    while parent[node] != -1:

        up = int(parent[node])

        length = from_a[node] - from_a[up]

        if from_a[up] <= half <= from_a[node]:

            return up, node, half - from_a[up]

        node = up

    # Degenerate case: put the root on the first edge.
    u, v, length = edges[0]

    return u, v, length / 2.0


def root_tree(edges, n_nodes, n_tips, mode, tip_names, outgroup):
    """
    Inserts a root node on one edge and returns
    (children, root, n_nodes) where children maps a node to a list of
    (child, branch length).
    """

    if mode == "none":

        # Leave the tree unrooted: hang it on the last internal node.
        neighbours = adjacency(edges, n_nodes)

        root = max(range(n_tips, n_nodes), key=lambda n: len(neighbours[n]))

        children = {}

        stack = [(root, -1)]

        while stack:

            node, came_from = stack.pop()

            children[node] = []

            for neighbour, length in neighbours[node]:

                if neighbour == came_from:
                    continue

                children[node].append((neighbour, length))

                stack.append((neighbour, node))

        return children, root, n_nodes


    if mode == "outgroup":

        if outgroup is None:
            raise ValueError("--root outgroup needs --outgroup")

        if outgroup not in tip_names:
            raise ValueError(
                f"Outgroup '{outgroup}' is not one of the samples"
            )

        tip = tip_names.index(outgroup)

        edge = next(
            (u, v, length)
            for u, v, length in edges
            if tip in (u, v)
        )

        u, v, length = edge

        # Root halfway along the outgroup's own branch.
        if u == tip:
            u, v = v, u

        anchor_u, anchor_v, offset = u, v, length / 2.0

    else:

        anchor_u, anchor_v, offset = midpoint_of_tree(
            edges, n_nodes, n_tips
        )


    root = n_nodes

    rebuilt = [
        (u, v, length)
        for u, v, length in edges
        if {u, v} != {anchor_u, anchor_v}
    ]

    total = next(
        length
        for u, v, length in edges
        if {u, v} == {anchor_u, anchor_v}
    )

    rebuilt.append((root, anchor_u, offset))
    rebuilt.append((root, anchor_v, max(total - offset, 0.0)))

    n_nodes += 1

    neighbours = adjacency(rebuilt, n_nodes)

    children = {}

    stack = [(root, -1)]

    while stack:

        node, came_from = stack.pop()

        children[node] = []

        for neighbour, length in neighbours[node]:

            if neighbour == came_from:
                continue

            children[node].append((neighbour, length))

            stack.append((neighbour, node))

    return children, root, n_nodes


# --------------------------------------------------------------------- #
#  Newick
# --------------------------------------------------------------------- #

def to_newick(children, root, n_tips, tip_names, support=None) -> str:

    def render(node):

        if node < n_tips:

            name = tip_names[node].replace(" ", "_")

            return name

        parts = []

        for child, length in children[node]:

            parts.append(f"{render(child)}:{length:.6f}")

        label = ""

        if support is not None and node in support:

            label = f"{support[node]:.0f}"

        return "(" + ",".join(parts) + ")" + label

    sys.setrecursionlimit(max(10_000, n_tips * 20))

    return render(root) + ";"


def support_per_node(children, root, n_tips, split_counts, replicates):
    """
    Maps each internal node to the bootstrap support of the split
    defined by the edge above it.
    """

    support = {}

    def descend(node):

        if node < n_tips:
            return {node}

        tips = set()

        for child, _ in children[node]:
            tips |= descend(child)

        if node != root:

            side = frozenset(tips)

            if 0 in side:
                side = frozenset(set(range(n_tips)) - side)

            if 1 < len(side) < n_tips - 1:

                support[node] = (
                    100.0 * split_counts.get(side, 0) / replicates
                )

        return tips

    descend(root)

    return support


# --------------------------------------------------------------------- #

def main():

    args = parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)


    genotypes, names = read_traw(args.traw)

    n_samples, n_loci = genotypes.shape

    print(f"{n_samples} samples x {n_loci} loci")


    g0, g1, g2, present = one_hot(genotypes)

    del genotypes

    distance = ibs_distance(g0, g1, g2, present)

    if np.isnan(distance).any():
        raise ValueError(
            "Some pairs share no called locus; filter the input further"
        )


    #
    # Cross-check against PLINK, if its matrix was supplied.
    #

    if args.mdist is not None and args.mdist.exists():

        plink = np.loadtxt(args.mdist)

        difference = np.abs(plink - distance).max()

        print(f"max |PLINK - computed here| = {difference:.6f}")

        if difference > 1e-4:
            raise ValueError(
                f"Distance matrix disagrees with PLINK by {difference:.6f}; "
                f"refusing to build a tree on it"
            )


    pd.DataFrame(
        distance,
        index=names,
        columns=names
    ).to_csv(
        args.outdir / "distance_matrix.tsv",
        sep="\t",
        float_format="%.6f"
    )


    edges, n_nodes = neighbour_joining(distance)

    print(f"NJ tree: {len(edges)} edges over {n_nodes} nodes")


    #
    # Bootstrap.
    #

    split_counts: dict[frozenset, int] = {}

    if args.bootstrap > 0:

        rng = np.random.default_rng(args.seed)

        for replicate in range(args.bootstrap):

            weights = rng.multinomial(
                n_loci,
                np.full(n_loci, 1.0 / n_loci)
            ).astype(np.float32)

            replicate_distance = ibs_distance(
                g0, g1, g2, present, weights
            )

            if np.isnan(replicate_distance).any():
                continue

            replicate_edges, replicate_nodes = neighbour_joining(
                replicate_distance
            )

            for side in splits_of(
                replicate_edges, replicate_nodes, n_samples
            ):
                split_counts[side] = split_counts.get(side, 0) + 1

            if (replicate + 1) % 10 == 0:
                print(f"  bootstrap {replicate + 1}/{args.bootstrap}")


    children, root, n_nodes = root_tree(
        edges,
        n_nodes,
        n_samples,
        args.root,
        names,
        args.outgroup
    )


    support = None

    if args.bootstrap > 0:

        support = support_per_node(
            children,
            root,
            n_samples,
            split_counts,
            args.bootstrap
        )

        values = np.array(list(support.values()))

        if len(values):
            print(
                f"bootstrap support: median {np.median(values):.0f}%, "
                f"{int((values >= 70).sum())} of {len(values)} "
                f"nodes at or above 70%"
            )


    newick = to_newick(children, root, n_samples, names, support)

    (args.outdir / "tree.nwk").write_text(newick + "\n")

    pd.DataFrame({"tip": names}).to_csv(
        args.outdir / "tree_tips.tsv",
        sep="\t",
        index=False
    )

    print()
    print(args.outdir / "tree.nwk")


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
