"""Refit retained GCM nuisance models and historical fitted contrasts."""
from pathlib import Path
import json
import os
import numpy as np
import pandas as pd
import fig5_estimators as model
import fig5_spec as spec

pkg = Path(__file__).resolve().parents[1]
output = Path(os.environ.get("NAKDONG_OUTPUT", str(pkg / "run"))).resolve()
target = output / "results/fig5_full_refits"
target.mkdir(parents=True, exist_ok=True)
model.OUTPUT = target
base = pkg / "data/derived"
ai = pd.read_csv(base / "optionC_ai_analysis_full.csv").set_index("sample_id")
manifest = pd.read_csv(base / "v10_site_aware_sensitivity/fig5_site_blocked/cohort_sample_manifest.csv")
cohorts = {name: ai.loc[m.sort_values("row_index_zero_based").sample_id].reset_index()
           for name, m in manifest.groupby("cohort", sort=False)}
model.run_gcm_site_blocked(spec, cohorts)
model.run_gcomp_site_blocked(spec, cohorts)
checks = {}
for filename, keys in [
    ("gcm_site_blocked.csv", ["cohort", "driver", "community_pcs_k"]),
    ("gcomputation_site_blocked.csv", ["cohort", "driver"])]:
    new = pd.read_csv(target / filename).set_index(keys).sort_index()
    old = pd.read_csv(base / "v10_site_aware_sensitivity/fig5_site_blocked" / filename).set_index(keys).sort_index()
    assert new.index.equals(old.index)
    errors = {}
    for col in new.select_dtypes(include="number").columns:
        if col in old:
            a = new[col].to_numpy(float); b = old[col].to_numpy(float)
            if not np.allclose(a, b, rtol=1e-8, atol=1e-8, equal_nan=True):
                errors[col] = float(np.nanmax(np.abs(a-b)))
    checks[filename] = {"rows": len(new), "numeric_mismatches": errors}
(output / "qc/full_fig5_refits.json").write_text(json.dumps(checks, indent=2))
print(json.dumps(checks, indent=2))

