#!/usr/bin/env python3
"""Inspect software versions; do not install dependencies or fit models."""
from __future__ import annotations
import argparse
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--rscript", default=shutil.which("Rscript"))
args = parser.parse_args()
packages = []
for line in (ROOT / "core/environment/requirements.txt").read_text().splitlines():
    if not line.strip() or line.startswith("#"):
        continue
    name, expected = line.split("==")
    try:
        actual = importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        actual = None
    packages.append({"package": name, "expected": expected, "installed": actual,
                     "matches_recorded_version": actual == expected})
r_result = {"available": False, "packages": [], "error": None}
if args.rscript:
    expression = '''cat(paste("R",as.character(getRversion()),sep="\\t"),"\\n",sep=""); for (p in c("DESeq2","vegan","permute","igraph")) {v <- if(requireNamespace(p,quietly=TRUE)) as.character(packageVersion(p)) else "MISSING";cat(paste(p,v,sep="\\t"),"\\n",sep="")}'''
    try:
        result = subprocess.run([args.rscript, "--vanilla", "-e", expression],
                                capture_output=True, text=True, timeout=60)
        r_result["available"] = result.returncode == 0
        r_result["packages"] = [dict(zip(["package", "installed"], line.split("\t", 1)))
                                for line in result.stdout.splitlines() if "\t" in line]
        if result.returncode:
            r_result["error"] = result.stderr[-2000:]
    except (OSError, subprocess.TimeoutExpired) as error:
        r_result["error"] = str(error)
expected_r = {"R": "4.4.3", "DESeq2": "1.46.0", "vegan": "2.7.2", "permute": "0.9.8"}
for record in r_result["packages"]:
    record["expected"] = expected_r.get(record["package"])
    record["matches_recorded_version"] = (record["installed"] == record["expected"]
                                           if record["expected"] else None)
ok = (all(p["matches_recorded_version"] for p in packages) and r_result["available"]
      and len(r_result["packages"]) == 5
      and all(p["installed"] != "MISSING" for p in r_result["packages"])
      and all(p["matches_recorded_version"] is not False for p in r_result["packages"]))
print(json.dumps({"status": "PASS" if ok else "REVIEW_ENVIRONMENT",
                  "python_version": platform.python_version(), "python_packages": packages,
                  "R": r_result, "new_model_fits": 0, "automatic_installs": 0}, indent=2))
raise SystemExit(0 if ok else 1)
