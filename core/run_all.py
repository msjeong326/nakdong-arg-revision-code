#!/usr/bin/env python3
"""Execute the portable core and write restartable status without external actions."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import traceback
import numpy as np
import pandas as pd

PKG = Path(__file__).resolve().parent
OUT = Path(os.environ.get("NAKDONG_OUTPUT", str(PKG/"run"))).resolve()
STAGES = ["01_arg_associations.py", "02_community_DESeq2.R", "03_functional_MAG.py",
          "04_DML_contrasts.py", "05_predictive_mirror.py", "06_spatial_network_procrustes.R"]

def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()

def write_state(data: dict) -> None:
    data["updated_utc"] = datetime.now(timezone.utc).isoformat()
    tmp=OUT/"status.json.tmp"; tmp.write_text(json.dumps(data,indent=2));tmp.replace(OUT/"status.json")

def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument("--rscript",default=shutil.which("Rscript"))
    parser.add_argument("--full-fig5",action="store_true")
    parser.add_argument("--resume",action="store_true")
    args=parser.parse_args()
    for name in ["results","qc","logs"]: (OUT/name).mkdir(parents=True,exist_ok=True)
    if not args.rscript or not Path(args.rscript).exists():
        raise SystemExit("Rscript is required; set --rscript to an R installation with the documented packages.")
    env=dict(os.environ,NAKDONG_OUTPUT=str(OUT),OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1")
    files=[p for d in ["data","code","expected"] for p in (PKG/d).rglob("*") if p.is_file() and "__pycache__" not in p.parts]
    before={str(p.relative_to(PKG)):sha(p) for p in files}
    fingerprint=hashlib.sha256(json.dumps(before,sort_keys=True).encode()).hexdigest()
    state={"state":"RUNNING","pid":os.getpid(),"package_fingerprint":fingerprint,
           "completed_stages":[],"failed_stages":[],"public_release":False}
    if args.resume and (OUT/"status.json").exists():
        old=json.loads((OUT/"status.json").read_text())
        if old.get("package_fingerprint")!=fingerprint:
            raise RuntimeError("Package changed; use a fresh output directory instead of --resume.")
        state["completed_stages"]=old.get("completed_stages",[])
    stages=STAGES+(["07_fig5_full_refits.py"] if args.full_fig5 else [])
    write_state(state)
    for stage in stages:
        if stage in state["completed_stages"]:continue
        state["active_stage"]=stage;write_state(state)
        command=[args.rscript if stage.endswith(".R") else sys.executable,str(PKG/"code"/stage)]
        print("START",stage,flush=True)
        with (OUT/"logs"/(stage+".log")).open("w") as log:
            result=subprocess.run(command,cwd=OUT,env=env,stdout=log,stderr=subprocess.STDOUT)
        if result.returncode:
            state["failed_stages"].append({"stage":stage,"exit_code":result.returncode})
            print("FAILED",stage,result.returncode,flush=True)
        else:
            state["completed_stages"].append(stage);print("PASS",stage,flush=True)
        write_state(state)
    comparisons=[]
    for ref in sorted((PKG/"expected").glob("*.csv")):
        got=OUT/"results"/ref.name
        record={"file":ref.name,"exists":got.exists()}
        if got.exists():
            a=pd.read_csv(ref);b=pd.read_csv(got)
            record["same_shape"]=a.shape==b.shape
            issues=[]
            if a.shape==b.shape and list(a.columns)==list(b.columns):
                for col in a:
                    if pd.api.types.is_numeric_dtype(a[col]) and pd.api.types.is_numeric_dtype(b[col]):
                        if not np.allclose(a[col],b[col],equal_nan=True,atol=1e-8,rtol=1e-8):issues.append(col)
                    elif not a[col].fillna("").equals(b[col].fillna("")):issues.append(col)
            else:issues.append("shape_or_columns")
            record["mismatching_columns"]=issues
            record["passed"]=not issues
        else:record["passed"]=False
        comparisons.append(record)
    preserved=all(sha(PKG/name)==digest for name,digest in before.items())
    validation={"comparisons":comparisons,"package_inputs_unchanged":preserved,
                "release_gate":"OPEN: factual input authority and complete manuscript coverage unresolved"}
    (OUT/"qc/package_validation.json").write_text(json.dumps(validation,indent=2))
    ok=not state["failed_stages"] and preserved and all(x["passed"] for x in comparisons)
    state.update(state="CORE_VALIDATED_RELEASE_PENDING" if ok else "VALIDATION_REQUIRES_REVIEW",
                 active_stage=None,validation=validation)
    write_state(state)
    print(json.dumps({"state":state["state"],"completed":len(state["completed_stages"]),
                      "failed":state["failed_stages"]},indent=2),flush=True)

if __name__=="__main__":
    try:main()
    except Exception:
        OUT.mkdir(parents=True,exist_ok=True)
        write_state({"state":"FAILED_RUNNER","error":traceback.format_exc(),"public_release":False})
        raise

