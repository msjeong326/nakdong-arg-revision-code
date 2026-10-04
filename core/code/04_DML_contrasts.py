# Portable adaptation of the audit implementation; statistical definitions unchanged.
from pathlib import Path
import os
PKG = Path(__file__).resolve().parents[1]
OUTPUT = Path(os.environ.get("NAKDONG_OUTPUT", str(PKG / "run"))).resolve()
for _folder in ["results", "qc"]:
    (OUTPUT / _folder).mkdir(parents=True, exist_ok=True)
from pathlib import Path
import ast,json
import numpy as np,pandas as pd
from scipy.stats import t
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold
from sklearn.metrics import r2_score
from statsmodels.stats.multitest import multipletests
from concurrent.futures import ThreadPoolExecutor
B=PKG/'data/derived';R=PKG/'data';O=OUTPUT
source=PKG/'code/fig5_spec.py'
tree=ast.parse(source.read_text());const={}
for node in tree.body:
 if isinstance(node,ast.Assign) and isinstance(node.targets[0],ast.Name):
  name=node.targets[0].id
  if name in ['EXPOSURES','CORE_CONTROLS','GCOMP_DRIVERS']:const[name]=ast.literal_eval(node.value)
pcs=[f'Comm_PC{i}' for i in range(1,11)]
features=pcs+['water_temp_C','MC_ppb','Chl_a_mg_m3','TOC_mg_L','NO3N_mg_L','discharge_mean_m3s','bloom','TN_mg_L','TP_mg_L','DO_mg_L_wq','pH_wq']
F=B/'v10_site_aware_sensitivity/fig5_site_blocked'
ai=pd.read_csv(B/'optionC_ai_analysis_full.csv').set_index('sample_id')
manifest=pd.read_csv(F/'cohort_sample_manifest.csv')
def run(item):
 cohort,man=item
 d=ai.loc[man.sort_values('row_index_zero_based').sample_id].copy();y=d.ARG_total.to_numpy(float);sites=d.site.to_numpy()
 rows=[]
 for treatment,(label,group) in const['EXPOSURES'].items():
  controls=[c for c in const['CORE_CONTROLS'] if c!=treatment]
  Z=np.column_stack([d[controls],pd.get_dummies(d.season,drop_first=True,dtype=float)])
  tr=d[treatment].to_numpy(float);tr=(tr-tr.mean())/tr.std()
  yr=np.zeros(len(d));trr=np.zeros(len(d))
  for train,test in GroupKFold(5,shuffle=True,random_state=42).split(Z,groups=sites):
   yp=RandomForestRegressor(n_estimators=120,random_state=42,n_jobs=1).fit(Z[train],y[train]).predict(Z[test])
   tp=RandomForestRegressor(n_estimators=120,random_state=42,n_jobs=1).fit(Z[train],tr[train]).predict(Z[test])
   yr[test]=y[test]-yp;trr[test]=tr[test]-tp
  denom=np.dot(trr,trr);theta=np.dot(trr,yr)/denom;score=trr*(yr-theta*trr)
  cs=np.array([score[sites==s].sum() for s in sorted(set(sites))])
  se=np.sqrt(np.dot(cs,cs))/denom;p=2*t.sf(abs(theta/se),len(cs)-1)
  rows.append({'cohort':cohort,'treatment':treatment,'theta':theta,'se':se,'p':p})
 dd=pd.DataFrame(rows);dd['q']=multipletests(dd.p,method='fdr_bh')[1]
 x=d[features].to_numpy(float)
 preds={col:np.zeros(len(d)) for _,col,_ in const['GCOMP_DRIVERS']}
 coherent=np.zeros(len(d));baseline=np.zeros(len(d));mismatches={}
 for train,test in GroupKFold(5,shuffle=True,random_state=7).split(x,groups=sites):
  rf=RandomForestRegressor(n_estimators=500,min_samples_leaf=3,random_state=0,n_jobs=1).fit(x[train],y[train])
  baseline[test]=rf.predict(x[test])
  for label,col,kind in const['GCOMP_DRIVERS']:
   j=features.index(col);lo,hi=(0.,1.) if kind=='binary' else np.quantile(d[col],[.1,.9])
   xl=x[test].copy();xh=xl.copy();xl[:,j]=lo;xh[:,j]=hi
   preds[col][test]=rf.predict(xh)-rf.predict(xl)
  jm=features.index('MC_ppb');jb=features.index('bloom');lo,hi=np.quantile(d.MC_ppb,[.1,.9])
  xl=x[test].copy();xh=xl.copy();xl[:,jm]=lo;xh[:,jm]=hi;xl[:,jb]=float(lo>1);xh[:,jb]=float(hi>1)
  coherent[test]=rf.predict(xh)-rf.predict(xl)
 gc=[{'cohort':cohort,'column':k,'shift_sd':v.mean()/y.std()} for k,v in preds.items()]
 gc.append({'cohort':cohort,'column':'MC_ppb_with_consistent_indicator_DIAGNOSTIC','shift_sd':coherent.mean()/y.std()})
 impossible=int(np.sum(d.bloom!=(d.MC_ppb>1)))
 proof={'cohort':cohort,'n':len(d),'raw_indicator_inconsistency':impossible,'binary_flip_counterfactuals_inconsistent':len(d),'binary_flip_total_counterfactuals':2*len(d),'p10_MC':float(np.quantile(d.MC_ppb,.1)),'p90_MC':float(np.quantile(d.MC_ppb,.9)),'crossfit_R2':r2_score(y,baseline)}
 print(cohort,'DONE',flush=True)
 return dd,pd.DataFrame(gc),proof
with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(run,list(manifest.groupby('cohort',sort=False))))
dml=pd.concat([v[0] for v in results]);gc=pd.concat([v[1] for v in results]);proof=[v[2] for v in results]
saved=pd.read_csv(F/'dml_site_blocked.csv');mm=dml.merge(saved,on=['cohort','treatment'],validate='one_to_one')
check={'DML_tests':len(mm),'theta_max_error':float(abs(mm.theta-mm.theta_z).max()),'se_max_error':float(abs(mm.se-mm.se_z_site_clustered).max()),'q_max_error':float(abs(mm.q-mm.q_bh24_t_cluster_df).max()),'gcomp_input_consistency':proof}
dml.to_csv(O/'results/independent_dml_refit.csv',index=False);gc.to_csv(O/'results/independent_gcomp_refit.csv',index=False)
(O/'qc/independent_ML_checks.json').write_text(json.dumps(check,ensure_ascii=False,indent=2))
print(json.dumps(check,ensure_ascii=False,indent=2));print(gc.to_string(index=False))

