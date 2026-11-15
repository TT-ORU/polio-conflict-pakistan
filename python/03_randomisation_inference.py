"""03  Randomisation inference for the interaction (20,000 permutations).

Each district's whole lagged-violence series is reassigned to another district
in the same province, the Mundlak terms are rebuilt from it, and the model is
refitted; the statistic is the cluster-robust t on the interaction, computed
identically in observed and permuted fits. About 7 minutes; FAST=1 uses 1,000.
"""
import numpy as np, pandas as pd
import statsmodels.api as sm, statsmodels.formula.api as smf
from common import (load, estimation_sample, prep, pois, Keys, TAB, KEYS,
                    F_INT, SEED, FAST)

K = Keys("03_randomisation")
rng = np.random.default_rng(SEED)
E = prep(estimation_sample(load()))
m = pois(F_INT, E)
obs_t = m.params["vm10:zd10"] / m.bse["vm10:zd10"]

ids = E.district_id.unique()
pid = E.province.values[np.unique(E.district_id.values, return_index=True)[1]]
V7 = E.conflict_viol_lag1.values.reshape(-1, 7).T          # 7 years x districts
B = 1000 if FAST else 20000
null_t = np.empty(B)
for b in range(B):
    o = np.arange(len(ids))
    for p in np.unique(pid):
        k = np.where(pid == p)[0]
        o[k] = rng.permutation(k)
    s = E.copy()
    s["conflict_viol_lag1"] = V7[:, o].T.ravel()
    s["vm"] = s.groupby("district_id").conflict_viol_lag1.transform("mean")
    s["vm10"] = s.vm / 10
    s["vd10"] = (s.conflict_viol_lag1 - s.vm) / 10
    mb = smf.glm(F_INT, data=s, family=sm.families.Poisson(),
                 offset=s.lpop).fit(cov_type="cluster",
                 cov_kwds={"groups": s.district_id})
    null_t[b] = mb.params["vm10:zd10"] / mb.bse["vm10:zd10"]

p = float(np.mean(np.abs(null_t) >= abs(obs_t)))
mc = np.sqrt(p * (1 - p) / B)
print(f"B = {B}; observed cluster-robust t = {obs_t:.3f}")
print(f"randomisation p = {p:.4f} (Monte Carlo SE {mc:.4f}; "
      f"95% MC interval {p-1.96*mc:.4f} to {p+1.96*mc:.4f})")
print("null |t| quantiles: 90% {:.2f}, 95% {:.2f}, 99% {:.2f}".format(
    *np.quantile(np.abs(null_t), [.9, .95, .99])))
K.add("ri_B", B)
K.add("ri_obs_t", obs_t)
K.add("ri_p", p)
K.add("ri_mc_se", mc)
np.save(f"{TAB}/ri_null_t.npy", null_t)
pd.DataFrame([dict(B=B, obs_t=obs_t, p=p, mc_se=mc)]).to_csv(
    f"{TAB}/randomisation_inference.csv", index=False)
K.save()
