#!/usr/bin/env python3
"""
Extract ADMIXTURE cross-validation errors from log files.

Input:
    log directory containing *.out

Output:
    TSV:

        K    CV

        2    0.52341
        3    0.50128
        ...

Author:
    Omics26 pipeline
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


CV_PATTERN = re.compile(
    r"CV error\s*\(K=(\d+)\)\s*:\s*([0-9.eE+-]+)"
)


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--logdir",
        required=True,
        type=Path,
        help="Directory containing ADMIXTURE log files"
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output TSV"
    )

    return parser.parse_args()


LOG_NAME_RE = re.compile(r"^log\d+\.out$")


def extract_cv(log_file: Path):

    text = log_file.read_text(errors="ignore")

    match = CV_PATTERN.search(text)

    if match is None:
        return None

    k = int(match.group(1))

    cv = float(match.group(2))

    return k, cv


def main():

    args = parse_args()

    if not args.logdir.exists():
        raise FileNotFoundError(args.logdir)

    rows = []

    # Glob broadly, then keep only canonical log<K>.out files
    files = sorted(args.logdir.glob("log*.out"))
    files = [f for f in files if LOG_NAME_RE.match(f.name)]

    if len(files) == 0:
        raise RuntimeError(
            f"No log<K>.out files found in {args.logdir}"
        )

    seen = set()
    skipped = 0

    for file in files:

        result = extract_cv(file)

        if result is None:
            skipped += 1
            print(
                f"  skip {file.name:20s}  (no CV line)",
                file=sys.stderr
            )
            continue

        k, cv = result

        if k in seen:
            raise RuntimeError(f"Duplicate K={k}")

        seen.add(k)

        rows.append(
            {
                "K": k,
                "CV": cv,
                "log_file": file.name
            }
        )

        print(
            f"{file.name:20s} "
            f"K={k:<2d} "
            f"CV={cv:.6f}"
        )

    if len(rows) == 0:
        raise RuntimeError(
            "No CV errors extracted from any log file"
        )

    df = (
        pd.DataFrame(rows)
        .sort_values("K")
        .reset_index(drop=True)
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(args.output, sep="\t", index=False)

    print()
    print(f"Written {len(df)} rows  (skipped {skipped} files)")
    print(args.output)


if __name__ == "__main__":

    try:
        main()

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)