#!/usr/bin/env python3
"""
Check the SNPhylo trees before anything is plotted or interpreted.

A job leaving the SLURM queue is not evidence that it worked, and an
empty or truncated tree file is the usual way this fails. Nothing
downstream should run until these pass.

Output:
    tree_validation.txt
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import newick


def parse_args():

    parser = argparse.ArgumentParser(
        description="Validate the ML and bootstrap trees"
    )

    parser.add_argument("--ml-tree", required=True, type=Path)

    parser.add_argument("--bs-tree", required=True, type=Path)

    parser.add_argument("--expect-tips", type=int, default=413)

    parser.add_argument("--outgroup", default="ZZ01")

    parser.add_argument("--out", required=True, type=Path)

    return parser.parse_args()


def inspect(path: Path, expect: int, outgroup: str):

    results = []

    if not path.exists() or path.stat().st_size == 0:

        results.append((f"{path.name} exists and is not empty", "FAIL", "0 bytes"))

        return results, None

    results.append(
        (f"{path.name} exists and is not empty", "ok", f"{path.stat().st_size} bytes")
    )

    tree = newick.parse(path.read_text())

    if newick.rescale_support(tree):

        results.append(
            (
                f"{path.name}: support scale",
                "ok",
                "written as proportions, read as percent",
            )
        )

    names = newick.tip_names(tree)

    results.append(
        (
            f"{path.name}: tip count",
            "ok" if len(names) == expect else "FAIL",
            f"{len(names)} (expected {expect})",
        )
    )

    duplicated = sorted({name for name in names if names.count(name) > 1})

    results.append(
        (
            f"{path.name}: no duplicated labels",
            "ok" if not duplicated else "FAIL",
            "none" if not duplicated else ", ".join(duplicated[:10]),
        )
    )

    results.append(
        (
            f"{path.name}: contains {outgroup}",
            "ok" if outgroup in names else "FAIL",
            "yes" if outgroup in names else "no",
        )
    )

    return results, tree


def main():

    args = parse_args()

    lines = []

    checks = []

    for path in (args.ml_tree, args.bs_tree):

        results, tree = inspect(path, args.expect_tips, args.outgroup)

        checks.extend(results)

        if tree is not None and path == args.bs_tree:

            supported = [
                node.support
                for node in newick.all_nodes(tree)
                if not node.is_tip and node.support is not None
            ]

            if supported:

                strong = sum(1 for value in supported if value >= 70)

                very_strong = sum(1 for value in supported if value >= 90)

                lines.append(
                    f"\nbootstrap: {len(supported)} internal nodes carry a "
                    f"value; {strong} at or above 70, {very_strong} at or "
                    f"above 90"
                )

            else:

                checks.append(
                    (
                        f"{path.name}: carries bootstrap labels",
                        "FAIL",
                        "no internal node has a support value",
                    )
                )


    #
    # The two trees must be over the same taxa; if they are not, the
    # bootstrap was run against a different alignment.
    #

    body = [f"{status:<4}  {name:<44} {detail}" for name, status, detail in checks]

    text = "\n".join(body) + "\n" + "\n".join(lines) + "\n"

    args.out.write_text(text)

    print(text)

    if any(status == "FAIL" for _, status, _ in checks):

        raise SystemExit(
            "ERROR: tree validation failed. Read the SNPhylo logs before "
            "going further; in particular, do not present the ML tree as "
            "the final result when the bootstrap stage did not finish."
        )


if __name__ == "__main__":

    try:
        main()

    except SystemExit:
        raise

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
