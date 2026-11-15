"""01  Case audit (Table S1), Table 2, the interaction and marginal effects
(Table 3, Table S3), and the rescaled IRRs behind the precision bounds."""
import numpy as np, pandas as pd
from scipy import stats
from common import (load, estimation_sample, prep, pois, ci, line, Keys, TAB,
                    F_MAIN, F_INT)

K = Keys("01_main_models")
d = load()

# ---- Table S1: case audit ----------------------------------------------------
yr = d.groupby("year").wpv1_cases.sum()
pv = d.groupby("province").wpv1_cases.sum().sort_values(ascending=False)
print("Table S1 case audit, 2018-2025")
print("  by year:    ", yr.to_dict())
print("  by province:", pv.to_dict(), "| total", int(yr.sum()))
K.add("audit_total", yr.sum())
for y, v in yr.items():
    K.add(f"audit_{y}", v)
pd.concat([yr.rename("cases").reset_index().assign(by="year"),
           pv.rename("cases").reset_index().assign(by="province")]
          ).to_csv(f"{TAB}/tableS1_case_audit.csv", index=False)

# ---- Table 2 -----------------------------------------------------------------
base = d[(d.year >= 2019) & d.conflict_viol_lag1.notna()
         & d.ghsl_population.notna()]
c1 = prep(base, drop_empty_prov=False)
c2 = prep(base)
c3 = prep(estimation_sample(d))
specs = [("c1", "Year FE only", "wpv1_cases ~ vm1 + vd1 + C(yr)", c1),
         ("c2", "Province + year FE", "wpv1_cases ~ vm1 + vd1 + C(prov) + C(yr)", c2),
         ("c3", "+ zero-dose", "wpv1_cases ~ vm1 + vd1 + zd10 + C(prov) + C(yr)", c3)]
rows = []
print("\nTable 2 (IRR per violent event; zero-dose per +10 pp)")
for k, lab, f, s in specs:
    m = pois(f, s)
    K.add(f"t2_{k}_n", len(s))
    print(f" ({k}) {lab}: {len(s)} district-years")
    for term, tl in [("vm1", "district mean"), ("vd1", "annual deviation"),
                     ("zd10", "zero-dose per 10 pp")]:
        if term in m.params:
            r = K.ci(f"t2_{k}_{term}", m, term)
            print(line(tl, *r))
            rows.append(dict(column=lab, term=tl, n=len(s), irr=r[0], lo=r[1],
                             hi=r[2], p=r[3]))
pd.DataFrame(rows).to_csv(f"{TAB}/table2.csv", index=False)

E = c3
K.add("n_districts", E.district_id.nunique())
K.add("n_cases", E.wpv1_cases.sum())
K.add("share_zero", (E.wpv1_cases == 0).mean())
K.add("n_dy_with_case", (E.wpv1_cases > 0).sum())

# ---- rescaled IRRs, precision bounds (column 2) -------------------------------
m2 = pois("wpv1_cases ~ vm10 + vd10 + C(prov) + C(yr)", c2)
print("\nRescaled, per +10 events/yr (column 2):")
for term, tl in [("vm10", "district mean"), ("vd10", "annual deviation")]:
    r = K.ci(f"per10_{term}", m2, term)
    print(line(tl, *r))
b_hi_geo = [v for k, v in K.rows if k == "per10_vm10_hi"][0]
b_hi_dyn = [v for k, v in K.rows if k == "per10_vd10_hi"][0]
print(f"  precision bounds: geographic <= +{(b_hi_geo-1)*100:.1f}%, "
      f"dynamic <= +{(b_hi_dyn-1)*100:.1f}% per 10 events/yr")

# ---- interaction and marginal effects (Table 3, Table S3) ----------------------
mi = pois(F_INT, E)
print("\nConflict x coverage interaction (per +10 events/yr x per +10 pp):")
print(line("interaction", *K.ci("int", mi, "vm10:zd10")))
V = mi.cov_params()
mrows = []
print("Marginal IRR of +10 violent events/yr, by zero-dose rate:")
for z in [5, 10, 20, 30, 40, 50, 60]:
    zz = z / 10
    b = mi.params["vm10"] + zz * mi.params["vm10:zd10"]
    se = np.sqrt(V.loc["vm10", "vm10"] + zz ** 2 * V.loc["vm10:zd10", "vm10:zd10"]
                 + 2 * zz * V.loc["vm10", "vm10:zd10"])
    r = (np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se),
         2 * stats.norm.sf(abs(b / se)))
    for s_, v in zip(["", "_lo", "_hi", "_p"], r):
        K.add(f"marg{z}{s_}", v)
    mrows.append(dict(zero_dose_pct=z, irr=r[0], lo=r[1], hi=r[2], p=r[3]))
    print(line(f"{z}%", *r))
pd.DataFrame(mrows).to_csv(f"{TAB}/table3_tableS3_marginal_effects.csv",
                           index=False)
K.save()
