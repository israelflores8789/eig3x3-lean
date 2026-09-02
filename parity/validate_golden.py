#!/usr/bin/env python3

# Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
# Released under Apache 2.0 license as described in the file LICENSE.
# Authors: Israel Flores-Arbolay

"""validate_golden.py — compare golden.json for drift

Used by the release workflow: the committed `generated/golden.json`
(expected) is compared against freshly regenerated output from
`just golden` (actual). Floats are compared with a combined
relative/absolute tolerance (numpy.allclose semantics); everything
else must match exactly. Exits 1 and prints mismatch paths on drift.

Meant to be used in the release workflow to catch golden file drift.
This script never modifies either file. Golden fixtures are updated
deliberately, locally, and reviewed in a PR -- never in automation.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from typing import Any

MAX_REPORT = 25


def compare(
    expected: Any,
    actual: Any,
    path: str,
    rtol: float,
    atol: float,
    mismatches: list[str],
) -> None:
    # bool is a subclass of int; guard it before the numeric branch so
    # that True vs 1 and True vs False are treated as structural changes.
    if isinstance(expected, bool) or isinstance(actual, bool):
        if expected != actual or type(expected) is not type(actual):
            mismatches.append(f"{path}: {expected!r} != {actual!r}")
        return

    if isinstance(expected, (int, float)) and isinstance(actual, (int, float)):
        e, a = float(expected), float(actual)
        if math.isnan(e) or math.isnan(a):
            if not (math.isnan(e) and math.isnan(a)):
                mismatches.append(f"{path}: {e} != {a}")
        elif math.isinf(e) or math.isinf(a):
            if e != a:
                mismatches.append(f"{path}: {e} != {a}")
        elif abs(e - a) > atol + rtol * abs(a):
            mismatches.append(f"{path}: {e} != {a} (|diff| = {abs(e - a):.3e})")
        return

    if isinstance(expected, dict) and isinstance(actual, dict):
        if expected.keys() != actual.keys():
            missing = sorted(expected.keys() - actual.keys())
            extra = sorted(actual.keys() - expected.keys())
            mismatches.append(
                f"{path}: key mismatch (missing={missing}, extra={extra})"
            )
            return
        for key in expected:
            compare(expected[key], actual[key], f"{path}.{key}", rtol, atol, mismatches)
        return

    if isinstance(expected, list) and isinstance(actual, list):
        if len(expected) != len(actual):
            mismatches.append(f"{path}: length {len(expected)} != {len(actual)}")
            return
        for i, (e_item, a_item) in enumerate(zip(expected, actual, strict=True)):
            compare(e_item, a_item, f"{path}[{i}]", rtol, atol, mismatches)
        return

    if expected != actual:
        mismatches.append(f"{path}: {expected!r} != {actual!r}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare two golden vector files with numeric tolerance."
    )
    parser.add_argument("expected", help="Committed golden.json")
    parser.add_argument("actual", help="Freshly regenerated golden.json")
    parser.add_argument(
        "--rtol", type=float, default=1e-12, help="Relative tolerance (default: 1e-12)"
    )
    parser.add_argument(
        "--atol", type=float, default=1e-14, help="Absolute tolerance (default: 1e-14)"
    )
    args = parser.parse_args()

    with open(args.expected) as f:
        expected = json.load(f)
    with open(args.actual) as f:
        actual = json.load(f)

    mismatches: list[str] = []
    compare(expected, actual, "$", args.rtol, args.atol, mismatches)

    if not mismatches:
        print(f"No drift: matches within rtol={args.rtol}, atol={args.atol}.")
        return 0

    print(
        f"Golden drift detected: {len(mismatches)} mismatch(es) "
        f"(rtol={args.rtol}, atol={args.atol}):",
        file=sys.stderr,
    )
    for line in mismatches[:MAX_REPORT]:
        print(f"  {line}", file=sys.stderr)
    if len(mismatches) > MAX_REPORT:
        print(f"  ... and {len(mismatches) - MAX_REPORT} more", file=sys.stderr)
    print(
        "Fix: regenerate deliberately with `just golden` in a PR and review the diff.",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
