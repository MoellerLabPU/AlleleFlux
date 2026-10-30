"""Compare any two AlleleFlux output trees with the regression contract.

The golden test compares a run against the committed snapshot.  This command
compares two runs you made yourself -- two branches, two commits, two configs
-- with the same comparator and the same FULL contract (``_scope.py``): every
tracked file, strict tolerance, the cache / profile / inputMetadata / figures
exclusions applied.

Usage::

    # run the pipeline on commit A and on commit B, each into its own folder
    python -m tests.regression.diff_runs \\
        /path/to/run_A/longitudinal /path/to/run_B/longitudinal

The first tree is the reference: a file present there but missing from the
second is reported as MISSING; an extra file in the second is tolerated.

Exit code 0 -> every reference file matched.  1 -> at least one differed; the
list is printed.  2 -> a path was not a directory.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from tests.regression._scope import FULL_EXCLUDE_SUBSTRINGS
from tests.regression.compare import compare_trees


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "reference", type=Path, help="Output tree taken as the truth (the golden side)."
    )
    parser.add_argument("candidate", type=Path, help="Output tree under test.")
    parser.add_argument(
        "--no-exclusions",
        action="store_true",
        help="Compare every tracked file, including the cache parquet, profile "
        "copies, inputMetadata and figures that the contract normally skips.",
    )
    args = parser.parse_args()

    for label, path in (("reference", args.reference), ("candidate", args.candidate)):
        if not path.is_dir():
            print(f"ERROR: {label} path is not a directory: {path}", file=sys.stderr)
            return 2

    fails = compare_trees(
        args.candidate,
        args.reference,
        exclude_substrings=None if args.no_exclusions else FULL_EXCLUDE_SUBSTRINGS,
    )
    if not fails:
        print("OK -- every reference file matched by the candidate.")
        return 0
    print(f"{len(fails)} file(s) differ:")
    for f in fails:
        print()
        print(f)
    return 1


if __name__ == "__main__":
    sys.exit(main())
