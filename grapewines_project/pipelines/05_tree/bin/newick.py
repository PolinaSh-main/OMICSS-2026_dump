#!/usr/bin/env python3
"""
Minimal Newick reader, re-rooter and writer.

Shared by validate_tree.py, annotate_tips.py and plot_tree.py so that
all three agree on what a tip is and where the root sits.

The one piece of real work here is re-rooting. A bootstrap value belongs
to a *split* of the tips, not to a node, but Newick can only write it as
a label on a node -- by convention, the node at the far end of the
branch that makes the split. Moving the root therefore has to move those
labels: on every branch whose direction is reversed, the value has to
travel to what is now the child. Re-rooting a tree by rewriting the
string, without doing this, silently shifts every support value by one
branch.
"""

from __future__ import annotations

import re
import sys


class Node:

    __slots__ = ("name", "length", "support", "children", "x", "y")

    def __init__(self, name="", length=0.0, support=None):

        self.name = name
        self.length = length
        self.support = support
        self.children = []

        self.x = 0.0
        self.y = 0.0

    @property
    def is_tip(self) -> bool:
        return not self.children


def parse(text: str) -> Node:
    """
    Recursive-descent parser. Accepts branch lengths and numeric
    internal labels, which are read as bootstrap support.
    """

    text = "".join(text.split())

    position = 0

    def parse_node() -> Node:

        nonlocal position

        node = Node()

        if position < len(text) and text[position] == "(":

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
                    f"Unexpected {text[position]!r} at offset {position}"
                )

        start = position

        while position < len(text) and text[position] not in "(),:;":
            position += 1

        label = text[start:position].strip().strip("'\"")

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

    sys.setrecursionlimit(200_000)

    return parse_node()


def all_nodes(root: Node) -> list[Node]:

    found = []

    stack = [root]

    while stack:

        node = stack.pop()

        found.append(node)

        stack.extend(node.children)

    return found


def tips(root: Node) -> list[Node]:
    """
    Tips in left-to-right drawing order.
    """

    found = []

    def walk(node):

        if node.is_tip:
            found.append(node)
            return

        for child in node.children:
            walk(child)

    sys.setrecursionlimit(200_000)

    walk(root)

    return found


def tip_names(root: Node) -> list[str]:

    return [tip.name for tip in tips(root)]


def rescale_support(root: Node) -> bool:
    """
    Put bootstrap support on the 0-100 scale.

    phangorn 2.12 writes proportions, older conventions and the "support
    >= 70" everyone quotes are percentages. A tree whose values are all
    at or below 1 is proportions -- read as percentages it would say
    every clade is unsupported, and the figure would simply come out
    blank with no error anywhere. bootstrap_tree.R already converts;
    this catches a tree that came from somewhere else.

    Returns True if anything was rescaled.
    """

    values = [
        node.support
        for node in all_nodes(root)
        if not node.is_tip and node.support is not None
    ]

    if not values or max(values) > 1.0:
        return False

    for node in all_nodes(root):

        if node.support is not None:
            node.support *= 100.0

    return True


# --------------------------------------------------------------------- #
#  Re-rooting
# --------------------------------------------------------------------- #

def _undirected(root: Node):
    """
    Adjacency of the underlying unrooted tree.

    A root of degree two is an artefact of the rooted representation --
    it is suppressed here, its two branches merged into one, so that
    re-rooting somewhere else does not leave a stray internal node with
    a single child.
    """

    adjacency: dict[int, list[tuple[Node, float, float | None]]] = {}

    nodes = all_nodes(root)

    for node in nodes:
        adjacency[id(node)] = []

    def link(a: Node, b: Node, length: float, support):

        adjacency[id(a)].append((b, length, support))
        adjacency[id(b)].append((a, length, support))

    for node in nodes:

        for child in node.children:

            if node is root and len(root.children) == 2:
                # Handled below, as a single merged branch.
                continue

            link(node, child, child.length, child.support)

    if len(root.children) == 2:

        left, right = root.children

        support = left.support if left.support is not None else right.support

        link(left, right, left.length + right.length, support)

        del adjacency[id(root)]

    return adjacency


def reroot_on_tip(root: Node, outgroup: str) -> Node:
    """
    Returns a new tree rooted on the branch leading to `outgroup`, with
    the root placed halfway along it.
    """

    target = None

    for tip in tips(root):

        if tip.name == outgroup:
            target = tip
            break

    if target is None:
        raise ValueError(f"'{outgroup}' is not a tip of this tree")

    adjacency = _undirected(root)

    neighbours = adjacency[id(target)]

    if len(neighbours) != 1:
        raise ValueError(
            f"'{outgroup}' has {len(neighbours)} neighbours; not a tip"
        )

    anchor, length, _ = neighbours[0]

    half = length / 2.0

    new_root = Node()

    def build(node: Node, came_from: Node, branch: float, support) -> Node:

        copy = Node(node.name, branch, support)

        for neighbour, edge_length, edge_support in adjacency[id(node)]:

            if neighbour is came_from:
                continue

            copy.children.append(
                build(neighbour, node, edge_length, edge_support)
            )

        return copy

    sys.setrecursionlimit(200_000)

    new_root.children.append(build(target, anchor, half, None))

    new_root.children.append(build(anchor, target, half, None))

    return new_root


# --------------------------------------------------------------------- #

def write(root: Node, decimals: int = 6) -> str:

    def render(node: Node) -> str:

        if node.is_tip:

            body = re.sub(r"[\s(),:;]", "_", node.name)

        else:

            body = (
                "("
                + ",".join(render(child) for child in node.children)
                + ")"
            )

            if node.support is not None:
                body += f"{node.support:g}"

        return f"{body}:{node.length:.{decimals}f}"

    sys.setrecursionlimit(200_000)

    inner = ",".join(render(child) for child in root.children)

    label = "" if root.support is None else f"{root.support:g}"

    return f"({inner}){label};"
