# Portable adaptation of the audit implementation; statistical definitions unchanged.
from pathlib import Path
import os
PKG = Path(__file__).resolve().parents[1]
OUTPUT = Path(os.environ.get("NAKDONG_OUTPUT", str(PKG / "run"))).resolve()
for _folder in ["results", "qc"]:
    (OUTPUT / _folder).mkdir(parents=True, exist_ok=True)
from pathlib import Path
import ast,types,json
import pandas as pd,numpy as np
B=PKG/'data/derived';R=PKG/'data';O=OUTPUT
target=O/'results/knockoff_rerun';target.mkdir(exist_ok=True)
path=PKG/'code/fig5_estimators.py';tree=ast.parse(path.read_text())
nodes=[n for n in tree.body if isinstance(n,(ast.Import,ast.ImportFrom,ast.FunctionDef))]
scope={'__name__':'audit_knockoff','__file__':str(path)}
exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),'exec'),scope)
scope.update(OUTPUT=target,N_SPLITS=5,KNOCKOFF_GROUP_SEED=17,KNOCKOFF_SEED_COUNT=30,KNOCKOFF_TREES=600)
old=ast.parse((PKG/'code/fig5_spec.py').read_text())
env=next(ast.literal_eval(n.value) for n in old.body if isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Name) and n.targets[0].id=='KNOCKOFF_ENV')
source=types.SimpleNamespace(KNOCKOFF_ENV=env,PCS=[f'Comm_PC{i}' for i in range(1,11)])
ai=pd.read_csv(B/'optionC_ai_analysis_full.csv').set_index('sample_id')
man=pd.read_csv(B/'v10_site_aware_sensitivity/fig5_site_blocked/cohort_sample_manifest.csv')
cohorts={c:ai.loc[x.sort_values('row_index_zero_based').sample_id].reset_index() for c,x in man.groupby('cohort',sort=False)}
scope['run_knockoff_site_blocked'](source,cohorts)
new=pd.read_csv(target/'knockoff_site_blocked_predictive_mirror.csv')
saved=pd.read_csv(B/'v10_site_aware_sensitivity/fig5_site_blocked/knockoff_site_blocked_predictive_mirror.csv')
m=new.merge(saved,on=['cohort','feature'],suffixes=('_new','_old'),validate='one_to_one')
cols=['heldout_W_mean','heldout_W_positive_stability','heldout_original_permutation_delta_mse_mean','heldout_knockoff_permutation_delta_mse_mean']
out={'rows':len(m),'method':'execution of source functions extracted by AST; no top-level original module execution; same source, not independent algorithm','max_error':{c:float(abs(m[c+'_new']-m[c+'_old']).max()) for c in cols},'top3':new.sort_values('heldout_W_mean',ascending=False).groupby('cohort').head(3)[['cohort','label','heldout_W_mean']].to_dict('records')}
(O/'qc/knockoff_reexecution.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))

