"""02  Robustness: Karachi sensitivity (Table S2), 2019-2023 and reservoir
restrictions, negative binomial, leave-one-province-out, leave-one-district-out,
stratification by zero-dose tertile (Table S4), and the tertile grid
(Figure 3B)."""
import numpy as np, pandas as pd
import statsmodels.api as sm
from common import (load, estimation_sample, prep, pois, ci, line, Keys, TAB,
                    F_INT)

K = Keys("02_robustness")
d = load()
base = d[(d.year >= 2019) & d.conflict_viol_lag1.notna()
         & d.ghsl_population.notna()]
F2 = "wpv1_cases ~ vm1 + vd1 + C(prov) + C(yr)"
rows = []

# ---- Table S2: Karachi, and other restrictions (column-2 specification) -------
kar = base.district_name.str.contains("Karachi", case=False, na=False)
agg = base.assign(district_id=np.where(kar, 9999, base.district_id))
agg = (agg.groupby(["district_id", "year", "province"], as_index=False)
       .agg(wpv1_cases=("wpv1_cases", "sum"),
            conflict_viol_lag1=("conflict_viol_lag1", "mean"),
            ghsl_population=("ghsl_population", "sum")))
variants = [("base", "Primary, province + year FE", base),
            ("nokar", "Excluding Karachi's four districts", base[~kar]),
            ("karcol", "Karachi collapsed into one unit", agg),
            ("y1923", "Restricted to 2019-2023", base[base.year <= 2023]),
            ("reserv", "Reservoir provinces only", base[base.province.isin(
                ["Khyber-Pakhtunkhwa", "Sindh", "Balochistan",
                 "Federally Administered Tribal Ar"])])]
print("Table S2 and restrictions (IRR per violent event, district mean / deviation)")
for k, lab, s in variants:
    s = prep(s)
    m = pois(F2, s)
    print(f" {lab}  (n = {len(s)})")
    for term, tl in [("vm1", "district mean"), ("vd1", "annual deviation")]:
        r = K.ci(f"rob_{k}_{term}", m, term)
        print(line(tl, *r))
        rows.append(dict(analysis=lab, term=tl, n=len(s), irr=r[0], lo=r[1],
                         hi=r[2], p=r[3]))

# negative binomial: alpha fixed at 1 as originally run, and alpha estimated
s = prep(base)
X = pd.get_dummies(s[["vm1", "vd1", "province", "year"]].astype(
    {"year": str}), columns=["province", "year"], drop_first=True, dtype=float)
X = sm.add_constant(X)
print(" Negative binomial (province + year FE)")
nb1 = sm.GLM(s.wpv1_cases, X, family=sm.families.NegativeBinomial(alpha=1.0),
             offset=s.lpop).fit(cov_type="cluster",
             cov_kwds={"groups": s.district_id})
nb2 = sm.NegativeBinomial(s.wpv1_cases, X, offset=s.lpop).fit(
    disp=0, maxiter=500, cov_type="cluster", cov_kwds={"groups": s.district_id})
for tag, m in [("nb_alpha1", nb1), ("nb_est", nb2)]:
    for term, tl in [("vm1", "district mean"), ("vd1", "annual deviation")]:
        r = K.ci(f"{tag}_{term}", m, term)
        lab = f"{'alpha = 1' if tag == 'nb_alpha1' else 'alpha estimated'}, {tl}"
        print(line(lab, *r))
        rows.append(dict(analysis=f"Negative binomial, {lab}", term=tl,
                         n=len(s), irr=r[0], lo=r[1], hi=r[2], p=r[3]))
print(f"  estimated alpha = {nb2.params['alpha']:.3f}")
K.add("nb_alpha_estimated", nb2.params["alpha"])
pd.DataFrame(rows).to_csv(f"{TAB}/tableS2_robustness.csv", index=False)

# ---- interaction robustness (Table S4) ----------------------------------------
e = estimation_sample(d)
E = prep(e)
srows = []
print("\nTable S4: interaction, leave one province out")
for p in sorted(E.province.unique()):
    s = prep(e[e.province != p])
    m = pois(F_INT, s)
    r = K.ci(f"lopo_{p[:8]}", m, "vm10:zd10")
    print(line(f"drop {p}", *r))
    srows.append(dict(analysis=f"Excluding {p}", irr=r[0], lo=r[1], hi=r[2], p=r[3]))

lodo = []
for i in E.district_id.unique():
    s = prep(e[e.district_id != i])
    m = pois(F_INT, s)
    lodo.append((i, *ci(m, "vm10:zd10")))
lodo = pd.DataFrame(lodo, columns=["district_id", "irr", "lo", "hi", "p"])
K.add("lodo_min", lodo.irr.min())
K.add("lodo_max", lodo.irr.max())
K.add("lodo_maxp", lodo.p.max())
print(f"  leave one district out, 120 refits: IRR {lodo.irr.min():.3f} to "
      f"{lodo.irr.max():.3f}, max p = {lodo.p.max():.4f}")
lodo.to_csv(f"{TAB}/tableS4_leave_one_district_out.csv", index=False)

dist = E.groupby("district_id").agg(
    vm=("vm", "first"), zd=("zero_dose_rate_mics", "first"),
    cases=("wpv1_cases", "sum"), py=("ghsl_population", "sum"),
    name=("district_name", "first"), province=("province", "first"))
dist["zt"] = pd.qcut(dist.zd, 3, labels=["low", "mid", "high"])
dist["vt"] = pd.qcut(dist.vm, 3, labels=["low", "mid", "high"])
print(" Stratified by zero-dose tertile, no product term:")
for t in ["low", "mid", "high"]:
    s = prep(e[e.district_id.isin(dist.index[dist.zt == t])])
    m = pois("wpv1_cases ~ vm10 + vd10 + C(prov) + C(yr)", s)
    r = K.ci(f"strat_{t}", m, "vm10")
    print(line(f"{t} zero-dose tertile ({s.district_id.nunique()} districts)", *r))
    srows.append(dict(analysis=f"Stratified, {t} zero-dose tertile",
                      irr=r[0], lo=r[1], hi=r[2], p=r[3]))
pd.DataFrame(srows).to_csv(f"{TAB}/tableS4_interaction_robustness.csv",
                           index=False)

# ---- tertile grid (Figure 3B) ------------------------------------------------
cells = []
for vt in ["low", "mid", "high"]:
    for zt in ["low", "mid", "high"]:
        s = dist[(dist.vt == vt) & (dist.zt == zt)]
        cells.append(dict(violence=vt, zero_dose=zt, districts=len(s),
                          cases=int(s.cases.sum()),
                          rate_per_million=s.cases.sum() / s.py.sum() * 1e6))
grid = pd.DataFrame(cells)
grid.to_csv(f"{TAB}/figure3B_tertile_grid.csv", index=False)
hh = dist[(dist.vt == "high") & (dist.zt == "high")]
rest = dist.drop(hh.index)
r_hh = hh.cases.sum() / hh.py.sum() * 1e6
r_rest = rest.cases.sum() / rest.py.sum() * 1e6
K.add("grid_hh_districts", len(hh))
K.add("grid_hh_cases", hh.cases.sum())
K.add("grid_hh_rate", r_hh)
K.add("grid_rest_rate", r_rest)
K.add("grid_fold", r_hh / r_rest)
K.add("grid_max_other", grid.rate_per_million[
    ~((grid.violence == "high") & (grid.zero_dose == "high"))].max())
print(f"\nFigure 3B: high/high cell {len(hh)} districts, {hh.cases.sum()} cases, "
      f"{r_hh:.2f} per million vs {r_rest:.2f} elsewhere ({r_hh/r_rest:.1f}-fold)")
hh[["name", "province", "cases"]].sort_values("province").to_csv(
    f"{TAB}/high_violence_high_zero_dose_districts.csv")
print("  by province:", hh.province.value_counts().to_dict())
K.save()
