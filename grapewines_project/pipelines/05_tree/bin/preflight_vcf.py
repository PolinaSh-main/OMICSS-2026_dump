#!/usr/bin/env python3
"""
Checks the rooted VCF before anything expensive is submitted.

The tree job runs for hours; every one of these checks costs seconds and
each of them corresponds to a way the run has already gone wrong for
somebody. A mismatch is a hard stop, not a warning: the task brief is
explicit that the input must not be swapped for a more convenient file.

Output:
    preflight.txt   the counts that were checked, for the run record
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path


def parse_args():

    parser = argparse.ArgumentParser(
        description="Validate the rooted tree input VCF"
    )

    parser.add_argument("--vcf", required=True, type=Path)

    parser.add_argument("--expect-samples", type=int, default=413)

    parser.add_argument("--expect-variants", type=int, default=119195)

    parser.add_argument("--outgroup", default="ZZ01")

    parser.add_argument("--out", required=True, type=Path)

    return parser.parse_args()


def open_vcf(path: Path):

    if path.suffix == ".gz":
        return gzip.open(path, "rt")

    return path.open()


def scan(path: Path) -> tuple[list[str], int, set[str]]:
    """
    One pass: sample names from the header, then count variant records
    and collect the chromosomes they sit on.
    """

    samples: list[str] = []

    variants = 0

    chromosomes: set[str] = set()

    with open_vcf(path) as handle:

        for line in handle:

            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue

            variants += 1

            chromosomes.add(line.split("\t", 1)[0])

    return samples, variants, chromosomes


def main():

    args = parse_args()

    if not args.vcf.exists():
        raise SystemExit(f"ERROR: no such file: {args.vcf}")

    samples, variants, chromosomes = scan(args.vcf)

    outgroup_count = samples.count(args.outgroup)

    duplicates = sorted(
        {name for name in samples if samples.count(name) > 1}
    )

    checks = [
        (
            "sample count",
            len(samples),
            args.expect_samples,
            len(samples) == args.expect_samples,
        ),
        (
            f"{args.outgroup} occurrences",
            outgroup_count,
            1,
            outgroup_count == 1,
        ),
        (
            "variant records",
            variants,
            args.expect_variants,
            variants == args.expect_variants,
        ),
        (
            "duplicate sample names",
            len(duplicates),
            0,
            not duplicates,
        ),
    ]

    lines = [
        f"input      {args.vcf}",
        f"ingroup    {len(samples) - outgroup_count} accessions",
        f"chromosomes {len(chromosomes)}: "
        + ",".join(sorted(chromosomes, key=lambda c: (len(c), c))),
        "",
    ]

    for name, got, expected, passed in checks:

        lines.append(
            f"{'ok  ' if passed else 'FAIL'}  {name:<24} "
            f"{got}  (expected {expected})"
        )

    if duplicates:
        lines.append(f"      duplicated: {', '.join(duplicates[:10])}")

    text = "\n".join(lines) + "\n"

    args.out.write_text(text)

    print(text)

    if not all(passed for _, _, _, passed in checks):

        raise SystemExit(
            "\nERROR: preflight failed. Do not substitute another VCF -- "
            "record the counts above and resolve the "
            "discrepancy."
        )


if __name__ == "__main__":

    try:
        main()

    except SystemExit:
        raise

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
