# Portable adaptation of the audit implementation; statistical definitions unchanged.
from pathlib import Path
import os
PKG = Path(__file__).resolve().parents[1]
OUTPUT = Path(os.environ.get("NAKDONG_OUTPUT", str(PKG / "run"))).resolve()
for _folder in ["results", "qc"]:
    (OUTPUT / _folder).mkdir(parents=True, exist_ok=True)
from pathlib import Path
import json, hashlib,itertools,platform
import numpy as np,pandas as pd,statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from scipy import stats
R=PKG/'data';B=PKG/'data/derived';O=OUTPUT
out={};inputs=[]
def read(p,**kw):
 p=Path(p);inputs.append({'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()});return pd.read_csv(p,**kw)
raw=read(R/'rgi_combined_results.csv')
raw_file_rows=len(raw)
raw=raw.loc[(raw.site!='UNKNOWN')&(raw.month!='UNKNOWN')&raw.Cut_Off.notna()].copy()
out['rgi_excluded_non_sample_rows']=raw_file_rows-len(raw)
old=read(R/'functions/arg_class_abundance.csv').set_index('sample_id')
repl=[("aminoglycoside antibiotic","aminoglycoside"),("tetracycline antibiotic","tetracycline"),("penicillin beta-lactam","beta_lactam"),("carbapenem","carbapenem"),("cephalosporin","cephalosporin"),("glycopeptide antibiotic","glycopeptide"),("macrolide antibiotic","macrolide"),("fluoroquinolone antibiotic","fluoroquinolone"),("rifamycin antibiotic","rifamycin"),("phenicol antibiotic","phenicol"),("sulfonamide antibiotic","sulfonamide"),("trimethoprim antibiotic","trimethoprim"),("diaminopyrimidine antibiotic","trimethoprim"),("peptide antibiotic","peptide"),("lincosamide antibiotic","lincosamide"),("phosphonic acid antibiotic","fosfomycin"),("fosfomycin","fosfomycin"),("oxazolidinone antibiotic","oxazolidinone"),("mupirocin antibiotic","mupirocin"),("streptogramin antibiotic","streptogramin")]
def cls(s):return next((v for k,v in repl if k in s),'other')
raw['origin']=raw.site+raw.month
raw['class_first']=raw.Drug_Class.fillna('').str.split(';').str[0].str.strip().map(cls)
recon=raw.groupby(['origin','class_first']).size().unstack(fill_value=0).reindex(index=old.index,columns=old.columns,fill_value=0)
out['rgi']={'rows':len(raw),'mags':raw.sample_id.nunique(),'origins':raw.origin.nunique(),'cutoffs':raw.Cut_Off.value_counts().to_dict(),'class_matrix_shape':list(old.shape),'class_matrix_max_error':int(abs(recon-old).to_numpy().max()),'multilabel_rows':int(raw.Drug_Class.fillna('').str.contains(';').sum())}
meta=read(B/'optionC_metadata_preceding30d.csv').set_index('sample_id')
ai=read(B/'optionC_ai_analysis_full.csv').set_index('sample_id')
out['arg_total']={'max_error_class_sum_vs_Fig5':int((old.sum(axis=1).reindex(ai.index)-ai.ARG_total).abs().max()),'zero_rows':int((ai.ARG_total==0).sum()),'n':len(ai)}
d=meta.join(old,how='inner')
X=pd.concat([pd.Series(1.,index=d.index,name='const'),pd.get_dummies(d.season,prefix='season',drop_first=True,dtype=float),((d.water_temp_C-d.water_temp_C.mean())/d.water_temp_C.std()).rename('temp'),(d.MC_ppb>1.).astype(float).rename('high'),pd.get_dummies(d.site,prefix='site',drop_first=True,dtype=float)],axis=1)
saved=read(B/'v10_site_aware_sensitivity/fig4_feature_screens/fig4c_broad_site_fixed_comparison.csv')
rows=[]
for feature in old.columns:
 y=np.log2(d[feature]+1.);fit=sm.OLS(y,X).fit();ci=fit.conf_int().loc['high']
 site=fit.get_robustcov_results(cov_type='cluster',groups=d.site,use_correction=True,df_correction=True,use_t=True)
 roundfit=fit.get_robustcov_results(cov_type='cluster',groups=d.month,use_correction=True,df_correction=True,use_t=True)
 j=list(X.columns).index('high')
 rows.append({'feature':feature,'beta':fit.params.high,'se':fit.bse.high,'p':fit.pvalues.high,'CI_low':ci.iloc[0],'CI_high':ci.iloc[1],'p_cluster_site':site.pvalues[j],'p_cluster_round':roundfit.pvalues[j]})
z=pd.DataFrame(rows);z['q18']=multipletests(z.p,method='fdr_bh')[1];z['q18_cluster_site']=multipletests(z.p_cluster_site,method='fdr_bh')[1];z['q18_cluster_round']=multipletests(z.p_cluster_round,method='fdr_bh')[1]
merge=z.merge(saved[saved.source=='ARG drug class'],on='feature',validate='one_to_one')
out['class_ols']={'n':len(d),'rank':int(np.linalg.matrix_rank(X)),'columns':X.shape[1],'high_n':int((d.MC_ppb>1).sum()),'beta_max_error':float(abs(merge.beta-merge.site_fixed_coefficient).max()),'p_max_error':float(abs(merge.p-merge.site_fixed_p_value).max()),'nominal_sig':int((z.p<.05).sum()),'q18_sig':int((z.q18<.05).sum()),'site_cluster_q18_sig':int((z.q18_cluster_site<.05).sum()),'round_cluster_q18_sig':int((z.q18_cluster_round<.05).sum()),'min_q18':z.q18.min()}
g=saved.copy();g['recomputed_q']=multipletests(g.site_fixed_p_value,method='fdr_bh')[1]
out['broad_BH']={'features':len(g),'q_max_error':float(abs(g.recomputed_q-g.site_fixed_bh_q_global_1064).max()),'significant':g[g.recomputed_q<.05][['source','feature','recomputed_q']].to_dict('records')}
z.to_csv(O/'results/independent_class_models.csv',index=False)
F=B/'v10_site_aware_sensitivity/fig5_site_blocked'
fold=read(F/'site_blocked_fold_assignments.csv')
bad=[];summary=[]
for k,x in fold.groupby(['cohort','method']):
 ids=x.sample_id.value_counts()
 overlaps=[bool(set(v.train_sites.split(';'))&set(v.test_sites.split(';'))) for v in x.itertuples()]
 summary.append({'cohort':k[0],'method':k[1],'n':len(x),'folds':x.fold_one_based.nunique(),'unique_sample':x.sample_id.nunique(),'each_sample_once':bool((ids==1).all()),'overlap_rows':int(sum(overlaps)),'site_assigned_multiple_folds':int((x.groupby('site').fold_one_based.nunique()>1).sum())})
out['folds']=summary
scores=read(F/'gcm_site_cluster_scores.csv');gc=read(F/'gcm_site_blocked.csv')
signs=np.array(list(itertools.product([-1,1],repeat=10)))
records=[]
for k,x in scores.groupby(['cohort','driver','community_pcs_k']):
 s=x.sort_values('site').site_residual_product_sum.to_numpy()
 ref=abs(s.sum());sim=abs(signs@s)
 # floating point tolerance only for exact equality; includes all-positive and all-negative.
 p=float(np.mean(sim>=ref-max(1.,ref)*1e-12))
 records.append({'cohort':k[0],'driver':k[1],'community_pcs_k':k[2],'recomputed_p':p})
zz=pd.DataFrame(records)
zz['recomputed_q']=zz.groupby(['cohort','community_pcs_k']).recomputed_p.transform(lambda x:multipletests(x,method='fdr_bh')[1])
mm=zz.merge(gc,on=['cohort','driver','community_pcs_k'],validate='one_to_one')
out['gcm']={'tests':len(mm),'score_rows':len(scores),'max_p_error':float(abs(mm.recomputed_p-mm.p_exact_site_signflip_2s).max()),'max_q_error':float(abs(mm.recomputed_q-mm.q_bh23_site_signflip).max()),'bloom_n107':mm[(mm.driver=='Bloom')&mm.cohort.str.contains('n107')&mm.community_pcs_k.isin([0,1,10])][['community_pcs_k','recomputed_p','recomputed_q']].to_dict('records')}
mm.to_csv(O/'results/independent_gcm_reconstruction.csv',index=False)
out['environment']={'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__}
(O/'qc/independent_numeric_checks.json').write_text(json.dumps(out,indent=2,ensure_ascii=False))
(O/'qc/independent_numeric_input_hashes.json').write_text(json.dumps(inputs,indent=2,ensure_ascii=False))
print(json.dumps(out,indent=2,ensure_ascii=False))
