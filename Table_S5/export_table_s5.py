#!/usr/bin/env python3
"""Reproduce the approved Table S5 numeric presentation, without fitting models.

This portable adaptation uses the field mapping and exponentiation rules of the
unchanged W11 build_table_s5.py. It intentionally excludes cluster intervals,
18-class BH recalculation, other input branches, and alternative model grids.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

BASE = Path(__file__).resolve().parent
FIELDS = [
    "feature", "beta", "se_ols", "ci_ols_low", "ci_ols_high",
    "count_plus_one_ratio", "ratio_ci_ols_low", "ratio_ci_ols_high",
    "p_ols", "site_fixed_bh_q_global_1064", "mde80_log2_bonf_18",
]


def build_rows() -> list[dict[str, str | float]]:
    with (BASE / "inputs/approved_ols_fields.csv").open(newline="") as stream:
        source = list(csv.DictReader(stream))
    if len(source) != 18 or len({row["feature"] for row in source}) != 18:
        raise ValueError("Expected exactly 18 unique approved classes")
    rows = []
    for original in source:
        row: dict[str, str | float] = {"feature": original["feature"]}
        row.update({key: float(value) for key, value in original.items()
                    if key != "feature"})
        row["count_plus_one_ratio"] = 2.0 ** row["beta"]
        row["ratio_ci_ols_low"] = 2.0 ** row["ci_ols_low"]
        row["ratio_ci_ols_high"] = 2.0 ** row["ci_ols_high"]
        rows.append(row)
    return rows


def verify(rows: list[dict[str, str | float]]) -> dict:
    with (BASE / "expected/Table_S5_approved_numeric.csv").open(newline="") as stream:
        expected = list(csv.DictReader(stream))
    if len(rows) != len(expected):
        raise AssertionError("Unexpected expected-table length")
    checked = 0
    maximum = 0.0
    for row, reference in zip(rows, expected):
        if row["feature"] != reference["feature"]:
            raise AssertionError("Expected class order differs")
        for key in FIELDS[1:]:
            actual, target = float(row[key]), float(reference[key])
            maximum = max(maximum, abs(actual - target))
            if not math.isclose(actual, target, rel_tol=1e-12, abs_tol=1e-12):
                raise AssertionError((row["feature"], key, actual, target))
            checked += 1
    return {"status": "PASS", "numeric_cells_checked": checked,
            "maximum_absolute_error": maximum, "new_model_fits": 0,
            "new_inference_or_power_calculations": 0,
            "scope": "Presentation of locked approved OLS results only"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=BASE / "run/Table_S5.csv")
    args = parser.parse_args()
    rows = build_rows()
    report = verify(rows)
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite: {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    args.output.with_suffix(".validation.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
