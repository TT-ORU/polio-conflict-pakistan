"""04  Hurdle model (Table S6) and the temporal-ordering test of the attack
exposures (Table S5)."""
import numpy as np, pandas as pd
import statsmodels.api as sm, statsmodels.formula.api as smf
from statsmodels.discrete.truncated_model import TruncatedLFNegativeBinomialP
from common import load, estimation_sample, prep, pois, line, Keys, TAB, DATA

K = Keys("04_hurdle_temporal")
d = load()
e = estimation_sample(d)
E = prep(e)
rows = []

# ---- hurdle -----------------------------------------------------------------
# Part 1: any case. Year enters linearly because the single national case in
# 2021 separates the year indicators in a binary model.
E["any_case"] = (E.wpv1_cases > 0).astype(int)
p1 = smf.logit("any_case ~ vm10 + vd10 + zd10 + C(prov) + year", data=E).fit(
    disp=0, cov_type="cluster", cov_kwds={"groups": E.district_id})
p1i = smf.logit("any_case ~ vm10*zd10 + vd10 + C(prov) + year", data=E).fit(
    disp=0, cov_type="cluster", cov_kwds={"groups": E.district_id})
print(f"Hurdle part 1: any case ({E.any_case.sum()} of {len(E)} district-years), OR")
for t, lab in [("vm10", "violence, district mean, per 10/yr"),
               ("vd10", "violence, annual deviation, per 10/yr"),
               ("zd10", "zero-dose, per 10 pp")]:
    r = K.ci(f"h1_{t}", p1, t)
    print(line(lab, *r))
    rows.append(dict(part="1 any case (OR)", term=lab, est=r[0], lo=r[1], hi=r[2], p=r[3]))
r = K.ci("h1_int", p1i, "vm10:zd10")
print(line("interaction", *r))
rows.append(dict(part="1 any case (OR)", term="interaction", est=r[0], lo=r[1], hi=r[2], p=r[3]))

# Part 2: count given at least one; province FE only (126 observations)
pos = prep(e[e.wpv1_cases > 0])
X = sm.add_constant(pd.get_dummies(pos[["vm10", "vd10", "zd10", "prov"]],
                                   columns=["prov"], drop_first=True).astype(float))
p2 = TruncatedLFNegativeBinomialP(pos.wpv1_cases.values, X.values, truncation=0,
                                  offset=pos.lpop.values).fit(disp=0, maxiter=500)
names = list(X.columns)
print(f"Hurdle part 2: count given >0 ({len(pos)} district-years), "
      "zero-truncated NB, IRR")
for t, lab in [("vm10", "violence, district mean, per 10/yr"),
               ("vd10", "violence, annual deviation, per 10/yr"),
               ("zd10", "zero-dose, per 10 pp")]:
    k = names.index(t)
    b, se, pv = p2.params[k], p2.bse[k], p2.pvalues[k]
    r = (np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se), pv)
    for s_, v in zip(["", "_lo", "_hi", "_p"], r):
        K.add(f"h2_{t}{s_}", v)
    print(line(lab, *r))
    rows.append(dict(part="2 count given >0 (IRR)", term=lab, est=r[0], lo=r[1], hi=r[2], p=r[3]))
pd.DataFrame(rows).to_csv(f"{TAB}/tableS6_hurdle.csv", index=False)

# ---- temporal ordering ------------------------------------------------------
# Split at 2021. ACLED ends in 2024, so the late exposure window is 2022-2024.
d["half"] = np.where(d.year <= 2021, "early", "late")
# national totals with each ACLED event counted once (the panel repeats
# Karachi's city-wide events in each of its five districts)
_a = pd.read_csv(f"{DATA}/acled_attacks_district_year.csv")
_x = pd.read_csv(f"{DATA}/acled_to_panel_district_crosswalk.csv")
_a = _a[_a.district.isin(set(_x.acled_district.str.split().str.join(" ")))]
tot = pd.Series({"early": _a[_a.year.between(2018, 2021)].n_polio.sum(),
                 "late": _a[_a.year.between(2022, 2024)].n_polio.sum()})
K.add("attacks_early", tot["early"])
K.add("attacks_late", tot["late"])
trows = []
print(f"\nTemporal ordering (polio-campaign attacks: {tot['early']:.0f} in "
      f"2018-2021, {tot['late']:.0f} in 2022-2024)")
for expo, scale, lab in [("polio_attacks", 1, "polio-campaign attacks"),
                         ("hworker_attacks", 1, "health-worker attacks"),
                         ("conflict_viol", 10, "violent events, per 10/yr")]:
    for oh in ["late", "early"]:
        eh = "early" if oh == "late" else "late"
        x = d[d.half == eh].groupby("district_id")[expo].mean() / scale
        s = d[d.half == oh].copy()
        s["x"] = s.district_id.map(x)
        s = s.dropna(subset=["x", "ghsl_population"])
        s["lpop"] = np.log(s.ghsl_population)
        kp, ky = s.groupby("province").wpv1_cases.sum(), s.groupby("year").wpv1_cases.sum()
        s = s[s.province.isin(kp[kp > 0].index) & s.year.isin(ky[ky > 0].index)].copy()
        s["prov"] = s.province.astype(str).astype("category")
        s["yr"] = s.year.astype(str).astype("category")
        m = pois("wpv1_cases ~ x + C(prov) + C(yr)", s)
        r = K.ci(f"order_{expo}_{eh}_to_{oh}", m, "x")
        direction = f"{eh} exposure -> {oh} cases"
        print(line(f"{lab}: {direction}", *r))
        trows.append(dict(exposure=lab, direction=direction, irr=r[0], lo=r[1], hi=r[2], p=r[3]))
pd.DataFrame(trows).to_csv(f"{TAB}/tableS5_temporal_order.csv", index=False)
K.save()
