"""05  Temporal validation of the high-risk rule (Table S7) and its value over a
district's own case history, as counts and as incidence (Table S8).

The rule uses only 2018-2021 violence and the baseline MICS zero-dose rate; the
outcome is cases reported in 2022-2025."""
import numpy as np, pandas as pd
from common import (load, estimation_sample, pois, line, Keys, TAB,
                    classify_temporal)

K = Keys("05_temporal_validation")
d = load()
full_ids = estimation_sample(d).district_id.unique()
b = classify_temporal(d, full_ids)
h = d[d.year <= 2021]
b["hist_cases"] = h.groupby("district_id").wpv1_cases.sum()
b["hist_any"] = (b.hist_cases > 0).astype(int)
b["log_hist"] = np.log1p(b.hist_cases)
b["log_hist_rate"] = np.log1p(b.hist_cases /
                              h.groupby("district_id").ghsl_population.mean() * 1e6)
K.add("tv_flagged", b.high_risk.sum())
b.to_csv(f"{TAB}/tableS7_classification_by_district.csv")

ev = d[(d.year >= 2022) & d.ghsl_population.notna() & d.district_id.isin(b.index)].copy()
for c in ["high_risk", "hist_any", "log_hist", "log_hist_rate"]:
    ev[c] = ev.district_id.map(b[c])
ev["lpop"] = np.log(ev.ghsl_population)
kp, ky = ev.groupby("province").wpv1_cases.sum(), ev.groupby("year").wpv1_cases.sum()
ev = ev[ev.province.isin(kp[kp > 0].index) & ev.year.isin(ky[ky > 0].index)].copy()
ev["prov"] = ev.province.astype(str).astype("category")
ev["yr"] = ev.year.astype(str).astype("category")

g = ev.groupby("high_risk").agg(cases=("wpv1_cases", "sum"), py=("ghsl_population", "sum"))
g["rate"] = g.cases / g.py * 1e6
for k, tag in [(1, "flag"), (0, "rest")]:
    K.add(f"tv_cases_{tag}", g.loc[k, "cases"])
    K.add(f"tv_py_{tag}", g.loc[k, "py"] / 1e6)
    K.add(f"tv_rate_{tag}", g.loc[k, "rate"])
K.add("tv_cases_total", g.cases.sum())
K.add("tv_crude_rr", g.loc[1, "rate"] / g.loc[0, "rate"])
print(f"Table S7: {b.high_risk.sum()} districts flagged from 2018-2021 data; in "
      f"2022-2025 they carried {g.loc[1,'cases']:.0f} of {g.cases.sum():.0f} cases "
      f"({g.loc[1,'rate']:.2f} vs {g.loc[0,'rate']:.2f} per million, crude RR "
      f"{g.loc[1,'rate']/g.loc[0,'rate']:.1f})")

# overlap with the full-period classification
from common import prep
E = prep(estimation_sample(d))
q = E.groupby("district_id").agg(vm=("vm", "first"), zd=("zero_dose_rate_mics", "first"))
q["hr"] = ((pd.qcut(q.vm, 3, labels=False) == 2) & (pd.qcut(q.zd, 3, labels=False) == 2)).astype(int)
ov = int(((b.high_risk == 1) & (q.hr.reindex(b.index) == 1)).sum())
K.add("tv_overlap", ov)
print(f"  {ov} of {b.high_risk.sum()} also in the full-period high/high set")

specs = [
    ("m_py", "flag, province + year FE", "wpv1_cases ~ high_risk + C(prov) + C(yr)", "high_risk"),
    ("m_y", "flag, year FE only", "wpv1_cases ~ high_risk + C(yr)", "high_risk"),
    ("m_none", "province + year only", "wpv1_cases ~ C(prov) + C(yr)", None),
    ("m_hcnt", "prior burden alone, log(1+cases)", "wpv1_cases ~ log_hist + C(prov) + C(yr)", "log_hist"),
    ("m_any_flag", "flag + any prior case", "wpv1_cases ~ high_risk + hist_any + C(prov) + C(yr)", "high_risk"),
    ("m_any_hist", "  any prior case", "wpv1_cases ~ high_risk + hist_any + C(prov) + C(yr)", "hist_any"),
    ("m_cnt_flag", "flag + log(1+prior cases)", "wpv1_cases ~ high_risk + log_hist + C(prov) + C(yr)", "high_risk"),
    ("m_cnt_hist", "  log(1+prior cases)", "wpv1_cases ~ high_risk + log_hist + C(prov) + C(yr)", "log_hist"),
    ("m_rate_flag", "flag + log(1+prior cases per million)", "wpv1_cases ~ high_risk + log_hist_rate + C(prov) + C(yr)", "high_risk"),
    ("m_rate_hist", "  log(1+prior cases per million)", "wpv1_cases ~ high_risk + log_hist_rate + C(prov) + C(yr)", "log_hist_rate"),
    ("m_both_flag", "flag + both prior terms", "wpv1_cases ~ high_risk + log_hist + log_hist_rate + C(prov) + C(yr)", "high_risk"),
]
rows = []
print("\nTable S8: incremental value over the district's own history (2022-2025 outcome)")
for key, lab, f, term in specs:
    m = pois(f, ev)
    K.add(f"{key}_aic", m.aic)
    if term:
        r = K.ci(key, m, term)
        print(line(lab, *r) + f"   AIC {m.aic:.1f}")
        rows.append(dict(model=lab, irr=r[0], lo=r[1], hi=r[2], p=r[3], aic=m.aic))
    else:
        print(f"  {lab:<56} AIC {m.aic:.1f}")
        rows.append(dict(model=lab, aic=m.aic))
pd.DataFrame(rows).to_csv(f"{TAB}/tableS8_incremental_value.csv", index=False)

# the no-prior-case subgroup: reported as indicative only
sub = ev[ev.hist_any == 0]
gs = sub.groupby("high_risk").agg(c=("wpv1_cases", "sum"), p=("ghsl_population", "sum"))
print(f"\nNo prior case in 2018-2021: {sub.district_id.nunique()} districts, "
      f"{sub.wpv1_cases.sum():.0f} later cases; flagged "
      f"{gs.loc[1,'c']:.0f} cases ({gs.loc[1,'c']/gs.loc[1,'p']*1e6:.2f}/m) vs "
      f"{gs.loc[0,'c']:.0f} ({gs.loc[0,'c']/gs.loc[0,'p']*1e6:.2f}/m). Too few "
      "cases to model: indicative only.")
K.add("noprior_districts", sub.district_id.nunique())
K.add("noprior_rate_flag", gs.loc[1, "c"] / gs.loc[1, "p"] * 1e6)
K.add("noprior_rate_rest", gs.loc[0, "c"] / gs.loc[0, "p"] * 1e6)
K.save()
