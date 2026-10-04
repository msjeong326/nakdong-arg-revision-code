"""Exact original MDE helper, provided for inspection; not run in W12.

MDE in the approved table = ncp_for_power(92, 0.05 / 18) * se_ols.
This is a fixed-design OLS planning calculation, not observed BH power.
"""
from functools import lru_cache
from scipy import stats, optimize

@lru_cache(maxsize=None)
def ncp_for_power(df, alpha):
    cutoff=stats.t.ppf(1-alpha/2,df)
    # Squared noncentral t has noncentral F(1, df, ncp**2).
    # Avoid SciPy 1.17 nct opposite-tail NaNs, and its ncf.sf bug at nc=0.
    def power(ncp):
        return stats.f.sf(cutoff**2,1,df) if ncp==0 else stats.ncf.sf(cutoff**2,1,df,ncp*ncp)
    return optimize.brentq(lambda ncp: power(ncp)-.8,0,30)
