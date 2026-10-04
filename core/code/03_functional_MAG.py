# Portable adaptation of the audit implementation; statistical definitions unchanged.
from pathlib import Path
import os
PKG = Path(__file__).resolve().parents[1]
OUTPUT = Path(os.environ.get("NAKDONG_OUTPUT", str(PKG / "run"))).resolve()
for _folder in ["results", "qc"]:
    (OUTPUT / _folder).mkdir(parents=True, exist_ok=True)
from pathlib import Path
import json
import numpy as np,pandas as pd,statsmodels.api as sm
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
R=PKG/'data';B=PKG/'data/derived';O=OUTPUT
meta=pd.read_csv(B/'optionC_metadata_preceding30d.csv').set_index('sample_id')
saved=pd.read_csv(B/'v10_site_aware_sensitivity/fig4_feature_screens/fig4c_broad_site_fixed_comparison.csv')
files={'ARG drug class':'arg_class_abundance.csv','KEGG module':'kegg_module_per_sample.csv','KEGG pathway':'kegg_pathway_abundance.csv','N-cycling gene':'ncycling_gene_abundance.csv'}
def design(d):
 return pd.concat([pd.Series(1.,index=d.index,name='const'),pd.get_dummies(d.season,prefix='season',drop_first=True,dtype=float),((d.water_temp_C-d.water_temp_C.mean())/d.water_temp_C.std()).rename('temp'),(d.MC_ppb>1.).astype(float).rename('high'),pd.get_dummies(d.site,prefix='site',drop_first=True,dtype=float)],axis=1)
rows=[]
for family,file in files.items():
 features=saved.loc[saved.source==family,'feature'].tolist()
 data=pd.read_csv(R/'functions'/file).set_index('sample_id')[features].join(meta[['site','season','water_temp_C','MC_ppb']],how='inner')
 x=design(data)
 for feature in features:
  yy=pd.to_numeric(data[feature],errors='coerce')
  if yy.nunique(dropna=True)<=1:
   beta=0.;p=1.
  else:
   yy=np.log2(yy.fillna(yy.median())+1);fit=sm.OLS(yy,x).fit()
   beta=float(fit.params.high);p=float(fit.pvalues.high)
  rows.append({'source':family,'feature':feature,'n':len(data),'beta':beta,'p':p})
z=pd.DataFrame(rows);z['q']=multipletests(z.p,method='fdr_bh')[1]
m=z.merge(saved,on=['source','feature'],validate='one_to_one')
out={'broad_independent_models':{'n_features':len(z),'beta_max_error':float(abs(m.beta-m.site_fixed_coefficient).max()),'p_max_error':float(abs(m.p-m.site_fixed_p_value).max()),'q_max_error':float(abs(m.q-m.site_fixed_bh_q_global_1064).max()),'nominal_p05':int((z.p<.05).sum()),'BH_q05':int((z.q<.05).sum()),'largest_y_using_nominal':float(-np.log10(z.p.min())),'largest_y_using_BH':float(-np.log10(z.q.min()))}}
z.to_csv(O/'results/independent_all_1064_models.csv',index=False)
mags=pd.read_csv(R/'mag_table.csv')
raw=pd.read_csv(R/'rgi_combined_results.csv');raw=raw[(raw.site!='UNKNOWN')&(raw.month!='UNKNOWN')&raw.Cut_Off.notna()].copy()
counts=raw.groupby('sample_id').size().reindex(mags.representative_MAG_ID).to_numpy()
strict=raw[raw.Cut_Off=='Strict'].groupby('sample_id').size().reindex(mags.representative_MAG_ID,fill_value=0).to_numpy()
out['MAG']={'n':len(mags),'table_total_max_error':int(abs(counts-mags.rgi_candidate_matches).max()),'table_strict_max_error':int(abs(strict-mags.rgi_strict_matches).max()),'rho_genome_broad':spearmanr(mags.genome_size_bp,counts).statistic,'rho_genome_strict':spearmanr(mags.genome_size_bp,strict).statistic,'strict_fraction_percent':100*sum(strict)/sum(counts)}
agg=mags.groupby('origin_sample').agg(bp=('genome_size_bp','sum'),count=('rgi_candidate_matches','sum'),strict=('rgi_strict_matches','sum'),n_mags=('representative_MAG_ID','size'))
agg['count_per_Mbp']=agg['count']/(agg.bp/1e6)
out['MAG']['origin_rho_total_bases_count']=spearmanr(agg.bp,agg['count']).statistic
out['MAG']['origin_rho_total_bases_normalized']=spearmanr(agg.bp,agg.count_per_Mbp).statistic
ai=pd.read_csv(B/'optionC_ai_analysis_full.csv').set_index('sample_id')
ag=agg.join(ai[['Comm_PC1']])
out['MAG']['origin_rho_PC1_count']=spearmanr(ag.Comm_PC1,ag['count']).statistic
out['MAG']['origin_rho_PC1_count_per_Mbp']=spearmanr(ag.Comm_PC1,ag.count_per_Mbp).statistic
# Matrix-transform sensitivity: retain every recorded class label rather than just first label.
lookup=[("aminoglycoside antibiotic","aminoglycoside"),("tetracycline antibiotic","tetracycline"),("penicillin beta-lactam","beta_lactam"),("carbapenem","carbapenem"),("cephalosporin","cephalosporin"),("glycopeptide antibiotic","glycopeptide"),("macrolide antibiotic","macrolide"),("fluoroquinolone antibiotic","fluoroquinolone"),("rifamycin antibiotic","rifamycin"),("phenicol antibiotic","phenicol"),("sulfonamide antibiotic","sulfonamide"),("trimethoprim antibiotic","trimethoprim"),("diaminopyrimidine antibiotic","trimethoprim"),("peptide antibiotic","peptide"),("lincosamide antibiotic","lincosamide"),("phosphonic acid antibiotic","fosfomycin"),("fosfomycin","fosfomycin"),("oxazolidinone antibiotic","oxazolidinone"),("mupirocin antibiotic","mupirocin"),("streptogramin antibiotic","streptogramin")]
def classes(v):
 return sorted(set(next((short for full,short in lookup if full in p.strip()),'other') for p in str(v).split(';')))
raw['origin']=raw.site+raw.month
classsets={s:classes(s) for s in raw.Drug_Class.fillna('').unique()}
exploded=raw.assign(mapped=raw.Drug_Class.fillna('').map(classsets)).explode('mapped')
xmat=exploded.groupby(['origin','mapped']).size().unstack(fill_value=0)
primary=pd.read_csv(R/'functions/arg_class_abundance.csv').set_index('sample_id')
xmat=xmat.reindex(index=primary.index,columns=primary.columns,fill_value=0)
d=xmat.join(meta[['site','season','water_temp_C','MC_ppb']],how='inner');x=design(d)
tests=[]
for col in primary.columns:
 fit=sm.OLS(np.log2(d[col]+1),x).fit()
 tests.append({'feature':col,'beta':fit.params.high,'p':fit.pvalues.high})
tests=pd.DataFrame(tests);tests['q18']=multipletests(tests.p,method='fdr_bh')[1]
tests.to_csv(O/'results/diagnostic_multilabel_classes.csv',index=False)
out['multilabel_sensitivity']={'originally_multilabel_hits':int(raw.Drug_Class.fillna('').str.contains(';').sum()),'matrix_changed_cells':int((xmat!=primary).to_numpy().sum()),'scope':'same 18 class columns; all recorded labels, exploratory diagnostic only','q18_significant':int((tests.q18<.05).sum()),'min_q18':tests.q18.min()}
(O/'qc/functional_MAG_checks.json').write_text(json.dumps(out,indent=2,ensure_ascii=False));print(json.dumps(out,indent=2,ensure_ascii=False))

