"""10  Checks requested in the third review (pre-specified in
review3/PRESPECIFICATION_review3_checks.md, sha256 38d61ca0...).

None of this changes the primary analysis. Writes outputs/tables/review3_*.csv.
"""
import os, sys
import numpy as np, pandas as pd, geopandas as gpd
import statsmodels.api as sm, statsmodels.formula.api as smf
from patsy import dmatrices
from common import (load, estimation_sample, prep, pois, ci, line, TAB, DATA,
                    F_INT, classify_temporal)

SEED = 20261019
d = load()
E = prep(estimation_sample(d))
m0 = pois(F_INT, E)
T = "vm10:zd10"
print(line("Primary interaction (Table 3)", *ci(m0, T)))
out = {}

# ---------------------------------------------------------------- B1 sample
print("\n=== B1. districts in and out of the interaction sample ===")
ids = set(E.district_id)
a = d[(d.year >= 2019)].copy()
a["in"] = a.district_id.isin(ids)
g = a.groupby("district_id").agg(province=("province", "first"), inn=("in", "first"),
                                 cases=("wpv1_cases", "sum"), pop=("ghsl_population", "mean"),
                                 viol=("conflict_viol_lag1", "mean"),
                                 zd=("zero_dose_rate_mics", "first"))
g["inc"] = g.cases / (g["pop"] * 7) * 1e6
rows = []
for k, s in g.groupby("inn"):
    rows.append(dict(group="in interaction sample" if k else "outside",
                     districts=len(s), cases_2019_2025=s.cases.sum(),
                     share_of_cases=s.cases.sum() / g.cases.sum(),
                     population_m=s["pop"].sum() / 1e6,
                     cases_per_million_py=s.cases.sum() / (s["pop"].sum() * 7) * 1e6,
                     districts_with_any_case=(s.cases > 0).sum(),
                     violent_events_per_yr_median=s.viol.median(),
                     violent_events_per_yr_mean=s.viol.mean()))
b1 = pd.DataFrame(rows); print(b1.round(3).to_string(index=False))
print(" excluded districts by province:", g[~g.inn].province.value_counts().to_dict())
print(" excluded districts with cases:", g[(~g.inn) & (g.cases > 0)][["province", "cases"]].assign(
    name=d.groupby("district_id").district_name.first()).to_dict("records"))
b1.to_csv(f"{TAB}/review3_B1_sample.csv", index=False)

# ---------------------------------------------------------------- B2 targeting
print("\n=== B2. targeting rules: 2018-2021 information, 2022-2025 cases ===")
full_ids = E.district_id.unique()
c = classify_temporal(d, full_ids)
h = d[d.year <= 2021]
c["prior_cases"] = h.groupby("district_id").wpv1_cases.sum()
c["prior_rate"] = c.prior_cases / h.groupby("district_id").ghsl_population.mean()
ev = d[(d.year >= 2022) & d.district_id.isin(c.index)]
c["later_cases"] = ev.groupby("district_id").wpv1_cases.sum()
c["later_py"] = ev.groupby("district_id").ghsl_population.sum()
n_joint = int(c.high_risk.sum())
top = lambda v: c.index.isin(c[v].sort_values(ascending=False).index[:n_joint])
rules = {
    "(a) top violence tertile": (c.vt == "high").values,
    "(b) top zero-dose tertile": (c.zt == "high").values,
    "(c) both: the paper's rule": (c.high_risk == 1).values,
    "(d) any case 2018-2021": (c.prior_cases > 0).values,
    f"(e) {n_joint} districts with most cases 2018-2021": top("prior_cases"),
    f"(f) {n_joint} districts with highest incidence 2018-2021": top("prior_rate"),
    "(g) (c) or (d)": ((c.high_risk == 1) | (c.prior_cases > 0)).values,
    f"(h) {n_joint} districts with highest violence 2018-2021": top("vm_train"),
    f"(i) {n_joint} districts with highest zero-dose": top("zd"),
}
rows = []
for lab, f in rules.items():
    fl, ot = c[f], c[~f]
    rf = fl.later_cases.sum() / fl.later_py.sum() * 1e6
    ro = ot.later_cases.sum() / ot.later_py.sum() * 1e6
    rows.append(dict(rule=lab, districts=int(f.sum()),
                     population_share=fl.later_py.sum() / c.later_py.sum(),
                     later_cases=int(fl.later_cases.sum()),
                     share_of_later_cases=fl.later_cases.sum() / c.later_cases.sum(),
                     rate_flagged=rf, rate_other=ro, rate_ratio=rf / ro if ro > 0 else np.nan))
b2 = pd.DataFrame(rows); print(b2.round(3).to_string(index=False))
# what each single-factor rule adds to the joint rule
for lab, f in [("joint but no prior case", (c.high_risk == 1) & (c.prior_cases == 0)),
               ("prior case but not joint", (c.high_risk == 0) & (c.prior_cases > 0))]:
    print(f"  {lab}: {int(f.sum())} districts, {int(c[f].later_cases.sum())} later cases")
b2.to_csv(f"{TAB}/review3_B2_targeting.csv", index=False)

# ---------------------------------------------------------------- B3 inference
print("\n=== B3a. wild cluster score bootstrap, null imposed (Kline and Santos 2012) ===")
y, X = dmatrices(F_INT, E, return_type="dataframe")
j = list(X.columns).index(T)
Xr = X.drop(columns=T)
fr = sm.GLM(y, Xr, family=sm.families.Poisson(), offset=E.lpop).fit()
mu = fr.predict(Xr, offset=E.lpop).values
u = X.values * (y.values.ravel() - mu)[:, None]                  # scores at restricted fit
H = (X.values * mu[:, None]).T @ X.values                          # information
rest = [k for k in range(X.shape[1]) if k != j]
proj = H[j, rest] @ np.linalg.inv(H[np.ix_(rest, rest)])
s_i = u[:, j] - u[:, rest] @ proj                                  # efficient score
cl = E.district_id.values
sg = pd.Series(s_i).groupby(cl).sum().values
G = len(sg)
stat = sg.sum() / np.sqrt((sg ** 2).sum())
rng = np.random.default_rng(SEED)
B = 9999
res = {}
for wname in ["Rademacher", "Webb"]:
    if wname == "Rademacher":
        w = rng.choice([-1.0, 1.0], size=(B, G))
    else:
        w = rng.choice(np.array([-np.sqrt(1.5), -1, -np.sqrt(.5), np.sqrt(.5), 1, np.sqrt(1.5)]), size=(B, G))
    num = w @ sg
    den = np.sqrt(((w * sg) ** 2).sum(1))
    tb = num / den
    res[wname] = (np.mean(np.abs(tb) >= abs(stat)))
    print(f"  {wname:<10} score statistic {stat:.3f}; bootstrap p = {res[wname]:.4f} (B = {B})")
# conventional score test for comparison (asymptotic)
from scipy import stats
print(f"  asymptotic cluster-robust score test p = {2*stats.norm.sf(abs(stat)):.4f}; "
      f"Wald p (Table 3) = {m0.pvalues[T]:.4f}")

print("\n=== B3b. excluding outbreak years ===")
rows = [dict(model="primary", irr=ci(m0, T)[0], lo=ci(m0, T)[1], hi=ci(m0, T)[2], p=ci(m0, T)[3])]
for lab, yrs in [("excluding 2024", [2024]), ("excluding 2024 and 2025", [2024, 2025]),
                 ("excluding 2019 and 2020", [2019, 2020])]:
    s = prep(estimation_sample(d)[lambda x: ~x.year.isin(yrs)].pipe(lambda x: x))
    # Mundlak means recomputed within the reduced window (prep does this)
    m = pois(F_INT, s)
    r = ci(m, T); print(line(f"{lab} (n = {len(s)}, cases = {s.wpv1_cases.sum():.0f})", *r))
    rows.append(dict(model=lab, n=len(s), cases=s.wpv1_cases.sum(), irr=r[0], lo=r[1], hi=r[2], p=r[3]))
b3 = pd.DataFrame(rows)
b3 = pd.concat([b3, pd.DataFrame([dict(model="wild cluster score bootstrap, Rademacher", p=res["Rademacher"]),
                                  dict(model="wild cluster score bootstrap, Webb", p=res["Webb"])])])
b3.to_csv(f"{TAB}/review3_B3_inference.csv", index=False)

# ---------------------------------------------------------------- B4 covariates
print("\n=== B4. district covariates ===")
cov = pd.read_csv(f"{DATA}/district_covariates_review3.csv")   # travel time 2015, OSM, VIIRS
cov["log_fac"] = np.log1p(cov.osm_health_facility_count / cov.district_area_km2 * 1000)
cov["ntl"] = cov.viirs_mean_2019_2024
cov["travel_h"] = cov.travel_time_to_city_mean_min / 60
cov["log_ntl"] = np.log(cov.ntl + 0.01)
Ec = E.merge(cov[["district_id", "travel_h", "log_fac", "log_ntl"]], on="district_id", how="left")
assert Ec[["travel_h", "log_fac", "log_ntl"]].notna().all().all()
print("  correlations with zero-dose and district-mean violence (district level):")
dl = Ec.groupby("district_id")[["zd10", "vm10", "travel_h", "log_fac", "log_ntl"]].first()
print(dl.corr().loc[["zd10", "vm10"], ["travel_h", "log_fac", "log_ntl"]].round(2).to_string())
rows = []
for lab, add in [("primary", ""), ("+ travel time", " + travel_h"), ("+ facility density", " + log_fac"),
                 ("+ night-time lights", " + log_ntl"),
                 ("+ all three", " + travel_h + log_fac + log_ntl")]:
    m = pois(F_INT + add, Ec)
    r = ci(m, T); rz = ci(m, "zd10")
    print(line(f"{lab}: interaction", *r) + f" | zero-dose main {rz[0]:.2f} ({rz[1]:.2f}-{rz[2]:.2f})")
    rows.append(dict(model=lab, irr=r[0], lo=r[1], hi=r[2], p=r[3], zd_irr=rz[0], zd_lo=rz[1], zd_hi=rz[2]))
b4 = pd.DataFrame(rows); b4.to_csv(f"{TAB}/review3_B4_covariates.csv", index=False)
print(f"  pre-specified criterion (interaction > 1.05 with all three): "
      f"{'met' if b4.irr.iloc[-1] > 1.05 else 'NOT met'}")

# ---------------------------------------------------------------- B5 Conley
print("\n=== B5. spatial inference ===")
gd = gpd.read_file(f"{DATA}/pakistan_districts_gadm_level3.gpkg").to_crs(32642)
cen = gd.set_index("district_id").geometry.centroid
xy = np.c_[cen.x.values, cen.y.values] / 1000
idx = {k: i for i, k in enumerate(cen.index)}
P = np.array([xy[idx[k]] for k in E.district_id])
Dm = np.sqrt(((P[:, None, :] - P[None, :, :]) ** 2).sum(-1))
y0, X0 = dmatrices(F_INT, E, return_type="dataframe")
mu0 = m0.predict(E, offset=E.lpop).values
S = X0.values * (y0.values.ravel() - mu0)[:, None]
Hb = np.linalg.inv((X0.values * mu0[:, None]).T @ X0.values)
jj = list(X0.columns).index(T)
b = m0.params[T]
rows = []
same = (E.district_id.values[:, None] == E.district_id.values[None, :])
for kern in ["uniform (pre-specified)", "Bartlett"]:
    for cut in [100, 200, 300]:
        if kern.startswith("uniform"):
            K_ = ((Dm <= cut) | same).astype(float)
        else:
            K_ = np.maximum(1 - Dm / cut, 0); K_[same] = 1.0
        V = Hb @ (S.T @ K_ @ S) @ Hb
        v = V[jj, jj]
        if v <= 0:
            print(f"  {kern:<24} {cut} km: variance not positive (kernel not positive semi-definite)")
            rows.append(dict(variance=f"{kern}, {cut} km", note="variance not positive")); continue
        se = np.sqrt(v); p = 2 * stats.norm.sf(abs(b / se))
        print(f"  {kern:<24} {cut} km: IRR {np.exp(b):.3f} ({np.exp(b-1.96*se):.3f}-{np.exp(b+1.96*se):.3f}) p = {p:.4f}")
        rows.append(dict(variance=f"{kern}, {cut} km", irr=np.exp(b), lo=np.exp(b - 1.96 * se),
                         hi=np.exp(b + 1.96 * se), p=p))
# coarser spatial clusters: administrative divisions, wild score bootstrap (Webb)
div = gd.set_index("district_id").division
dv = E.district_id.map(div).values
sgd = pd.Series(s_i).groupby(dv).sum().values
Gd = len(sgd)
stat_d = sgd.sum() / np.sqrt((sgd ** 2).sum())
w = rng.choice(np.array([-np.sqrt(1.5), -1, -np.sqrt(.5), np.sqrt(.5), 1, np.sqrt(1.5)]), size=(B, Gd))
tb = (w @ sgd) / np.sqrt(((w * sgd) ** 2).sum(1))
p_div = np.mean(np.abs(tb) >= abs(stat_d))
md = smf.glm(F_INT, data=E.assign(dv=dv), family=sm.families.Poisson(), offset=E.lpop).fit(
    cov_type="cluster", cov_kwds={"groups": pd.factorize(dv)[0]})
r = ci(md, T)
print(line(f"clustered by division ({Gd} clusters), Wald", *r))
print(f"  clustered by division, wild score bootstrap (Webb): p = {p_div:.4f}")
rows += [dict(variance=f"division clusters ({Gd}), Wald", irr=r[0], lo=r[1], hi=r[2], p=r[3]),
         dict(variance=f"division clusters ({Gd}), wild score bootstrap (Webb)", p=p_div)]
pd.DataFrame(rows).to_csv(f"{TAB}/review3_B5_spatial.csv", index=False)

# ---------------------------------------------------------------- A1 zero-dose definition
print("\n=== A1. zero-dose measure: standard age group ===")
z = pd.read_csv(f"{DATA}/zero_dose_check_review4.csv").set_index("district_id")
rows = []
for lab, col in [("corrected, 12-23 months (primary)", "zdw_1223"), ("corrected, 12-35 months", "zdw_1235"),
                 ("corrected, 0-35 months", "zdw_0035"), ("previous definition, 0-35 months", "zdw_previous")]:
    s = E.copy()
    s["zd10"] = s.district_id.map(z[col]) * 10
    m = pois(F_INT, s); r = ci(m, T)
    mm = pois("wpv1_cases ~ vm10 + vd10 + zd10 + C(prov) + C(yr)", s); rz = ci(mm, "zd10")
    print(line(f"{lab}: interaction", *r) + f" | zero-dose (main model) {rz[0]:.2f} ({rz[1]:.2f}-{rz[2]:.2f})")
    rows.append(dict(definition=lab, irr=r[0], lo=r[1], hi=r[2], p=r[3], zd_irr=rz[0], zd_lo=rz[1], zd_hi=rz[2]))
pd.DataFrame(rows).to_csv(f"{TAB}/review3_A1_zero_dose_definition.csv", index=False)

# ---------------------------------------------------------------- A5 Moran 2021
from libpysal.weights import Queen
from esda.moran import Moran
gq = gpd.read_file(f"{DATA}/pakistan_districts_gadm_level3.gpkg").set_index("district_id")
g21 = E[E.year == 2021]
mu21 = m0.predict(g21, offset=g21.lpop)
r21 = (g21.wpv1_cases - mu21) / np.sqrt(mu21.clip(lower=1e-9))
wq = Queen.from_dataframe(gq.loc[g21.district_id].reset_index(), use_index=False); wq.transform = "r"
np.random.seed(20261016)
mi21 = Moran(r21.values, wq, permutations=999)
print(f"\n=== A5. Moran 2021: I = {mi21.I:.4f}, permutation SD {mi21.seI_sim:.4f}, "
      f"p_sim = {mi21.p_sim:.3f}, normal p = {mi21.p_norm:.2f}; cases = {g21.wpv1_cases.sum():.0f}")

# ---------------------------------------------------------------- keys for 99_verify
from common import Keys
K = Keys("10_review3")
def kci(prefix, r):
    for s_, v in zip(["", "_lo", "_hi", "_p"], [r.irr, r.lo, r.hi, r.p]):
        if pd.notna(v):
            K.add(prefix + s_, v)
K.add("moran21_p_norm", mi21.p_norm)
for r in pd.read_csv(f"{TAB}/review3_A1_zero_dose_definition.csv").itertuples():
    kci({"corrected, 12-23 months (primary)": "zd_primary", "corrected, 12-35 months": "zd_1235",
         "corrected, 0-35 months": "zd_0035", "previous definition, 0-35 months": "zd_previous"}[r.definition], r)
for r in pd.read_csv(f"{TAB}/review3_B3_inference.csv").itertuples():
    key = {"primary": "inf_primary", "excluding 2024": "inf_ex2024",
           "excluding 2024 and 2025": "inf_ex202425", "excluding 2019 and 2020": "inf_ex201920",
           "wild cluster score bootstrap, Rademacher": "inf_wb_rad",
           "wild cluster score bootstrap, Webb": "inf_wb_webb"}[r.model]
    kci(key, r)
for r in pd.read_csv(f"{TAB}/review3_B4_covariates.csv").itertuples():
    kci({"primary": "cov_none", "+ travel time": "cov_travel", "+ facility density": "cov_fac",
         "+ night-time lights": "cov_ntl", "+ all three": "cov_all"}[r.model], r)
for r in pd.read_csv(f"{TAB}/review3_B5_spatial.csv").itertuples():
    kci("sp_" + "".join(ch for ch in r.variance.lower() if ch.isalnum()), r)
for i, r in enumerate(pd.read_csv(f"{TAB}/review3_B2_targeting.csv").itertuples()):
    K.add(f"tg{i}_districts", r.districts); K.add(f"tg{i}_share", 100 * r.share_of_later_cases)
    K.add(f"tg{i}_popshare", 100 * r.population_share); K.add(f"tg{i}_rr", r.rate_ratio)
    K.add(f"tg{i}_cases", r.later_cases)
for r in pd.read_csv(f"{TAB}/review3_B1_sample.csv").itertuples():
    t_ = "in" if r.group.startswith("in") else "out"
    K.add(f"smp_{t_}_districts", r.districts); K.add(f"smp_{t_}_cases", r.cases_2019_2025)
    K.add(f"smp_{t_}_rate", r.cases_per_million_py); K.add(f"smp_{t_}_pop", r.population_m)
    K.add(f"smp_{t_}_viol", r.violent_events_per_yr_mean)
K.save()
