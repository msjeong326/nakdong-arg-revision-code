#!/usr/bin/env python3
"""Verify snapshot hashes and Python syntax without executing scientific modules."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
manifest = json.loads((ROOT / "SHA256SUMS.json").read_text())
failures = []
for relative, record in manifest.items():
    path = ROOT / relative
    if not path.is_file():
        failures.append({"file": relative, "problem": "missing"})
        continue
    if hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
        failures.append({"file": relative, "problem": "hash mismatch"})
python_sources = []
for path in sorted(ROOT.rglob("*.py")):
    if any(part in {".venv", "run", "results_local", "__pycache__"} for part in path.parts):
        continue
    relative = str(path.relative_to(ROOT))
    try:
        compile(path.read_bytes(), relative, "exec")
        python_sources.append(relative)
    except SyntaxError as error:
        failures.append({"file": relative, "problem": str(error)})
report = {"status": "PASS" if not failures else "FAIL", "hashed_files": len(manifest),
          "python_sources_compiled_without_execution": len(python_sources),
          "new_model_fits": 0, "failures": failures}
print(json.dumps(report, indent=2))
raise SystemExit(bool(failures))
