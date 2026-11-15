"""06  Predicted incidence at tertile medians, with intervals and contrasts
(Table 4).

Intervals: parametric bootstrap, 10,000 draws from the cluster-robust
covariance. A nonparametric district-cluster bootstrap fails to converge in
about 40% of resamples (a resample that loses a province's case-carrying
districts separates that province effect), and intervals from the survivors
would be selected on convergence, so it is not used."""
import numpy as np, pandas as pd, patsy
import statsmodels.api as sm
from common import load, estimation_sample, prep, Keys, TAB, F_INT, SEED

K = Keys("06_predicted_rates")
rng = np.random.default_rng(SEED)
E = prep(estimation_sample(load()))
y, X = patsy.dmatrices(F_INT, E, return_type="dataframe")
m = sm.GLM(y, X, family=sm.families.Poisson(), offset=E.lpop.values).fit(
    cov_type="cluster", cov_kwds={"groups": E.district_id})

q = E.groupby("district_id").agg(vm10=("vm10", "first"), zd10=("zd10", "first"))
vmed = q.groupby(pd.qcut(q.vm10, 3, labels=["low", "mid", "high"]), observed=True).vm10.median()
zmed = q.groupby(pd.qcut(q.zd10, 3, labels=["low", "mid", "high"]), observed=True).zd10.median()
SC = [("ll", "Low violence, low zero-dose", vmed["low"], zmed["low"]),
      ("hl", "High violence, low zero-dose", vmed["high"], zmed["low"]),
      ("lh", "Low violence, high zero-dose", vmed["low"], zmed["high"]),
      ("hh", "High violence, high zero-dose", vmed["high"], zmed["high"])]
print(f"Tertile medians: violence {vmed['low']*10:.1f} / {vmed['high']*10:.1f} "
      f"events/yr; zero-dose {zmed['low']*10:.1f}% / {zmed['high']*10:.1f}%")

pop = E.ghsl_population.values
draws = rng.multivariate_normal(m.params.values, np.asarray(m.cov_params()), size=10000)
tot = {}
rows = []
for key, lab, vv, zz in SC:
    s = E.copy()
    s["vm10"], s["zd10"], s["vd10"] = vv, zz, 0.0
    Xs = np.asarray(patsy.build_design_matrices([X.design_info], s,
                                                return_type="dataframe")[0])
    tot[key] = (np.exp(Xs @ draws.T) * pop[:, None]).sum(0)
    pt = (np.exp(Xs @ m.params.values) * pop).sum() / pop.sum() * 1e6
    lo, hi = np.percentile(tot[key] / pop.sum() * 1e6, [2.5, 97.5])
    for s_, v in zip(["", "_lo", "_hi"], [pt, lo, hi]):
        K.add(f"pred_{key}{s_}", v)
    rows.append(dict(scenario=lab, rate_per_million=pt, lo=lo, hi=hi))
    print(f"  {lab:<32} {pt:.2f} ({lo:.2f} to {hi:.2f}) per million person-years")
pd.DataFrame(rows).to_csv(f"{TAB}/table4_predicted_rates.csv", index=False)

print("Contrasts:")
for key, lab, a, b in [("c_viol_lowzd", "high vs low violence, low zero-dose", "hl", "ll"),
                       ("c_viol_highzd", "high vs low violence, high zero-dose", "hh", "lh"),
                       ("c_both", "both high vs both low", "hh", "ll")]:
    r = tot[a] / tot[b]
    med, lo, hi = np.median(r), *np.percentile(r, [2.5, 97.5])
    for s_, v in zip(["", "_lo", "_hi"], [med, lo, hi]):
        K.add(f"{key}{s_}", v)
    print(f"  {lab:<40} {med:.2f} ({lo:.2f} to {hi:.2f})")
K.save()
