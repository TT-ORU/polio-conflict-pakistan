"""11  Audit of the ACLED attack exposures (review 3, check A2).

The panel's polio-campaign and health-worker counts were keyword matches over all
ACLED event types, so they include protests (for example vaccinators demanding pay).
This script rebuilds both counts, restricts them to violent event types, and reruns
the temporal-ordering test of Table S5 with the restricted counts.
"""
import numpy as np, pandas as pd
from common import load, pois, line, TAB, DATA

import sys
# raw ACLED events are not distributed (ACLED terms); pass the file path
ev = pd.read_csv(sys.argv[1], low_memory=False)
N = ev.notes.fillna("").str.lower()
A = (ev.actor1.fillna("") + "|" + ev.assoc_actor_1.fillna("") + "|" + ev.actor2.fillna("")
     + "|" + ev.assoc_actor_2.fillna("")).str.lower()
ev["polio"] = N.str.contains("polio|vaccinat|immuni")                      # reproduces n_polio exactly
ev["hw"] = A.str.contains("health workers") | N.str.contains(
    r"health worker|lady health|vaccinat|doctor|nurse|paramedic")         # approximates n_hworker
VIOL = ["Battles", "Violence against civilians", "Explosions/Remote violence"]
ev["violent"] = ev.event_type.isin(VIOL)
ev["violent_or_riot"] = ev.event_type.isin(VIOL + ["Riots"])

# composition, 2018-2024
s = ev[ev.year.between(2018, 2024)]
comp = pd.concat({k: s[s[k]].event_type.value_counts() for k in ["polio", "hw"]}, axis=1).fillna(0).astype(int)
comp.loc["Total"] = comp.sum()
print("Event types of the flagged events, 2018-2024:\n", comp.to_string())
comp.to_csv(f"{TAB}/review3_A2_attack_event_types.csv")

# reproduce the file's counts
pol = pd.read_csv(f"{DATA}/acled_polio_district_year.csv")   # the earlier all-event counts
for k, col in [("polio", "n_polio"), ("hw", "n_hworker")]:
    mine = ev[ev[k]].groupby("year").size().loc[2015:2024]
    theirs = pol.groupby("year")[col].sum().loc[2015:2024]
    print(f"{col}: rebuilt/used ratio by year {np.round((mine / theirs).values, 2)}")

# district-year counts on the panel geography
ax = pd.read_csv(f"{DATA}/acled_to_panel_district_crosswalk.csv")
ax["acled_district"] = ax.acled_district.str.split().str.join(" ")
ev["district"] = ev.admin2.fillna("").str.upper().str.split().str.join(" ")
x = ev.merge(ax, left_on="district", right_on="acled_district")
cnt = {}
for k in ["polio", "hw"]:
    for scope in ["all", "violent", "violent_or_riot"]:
        m = x[k] if scope == "all" else (x[k] & x[scope])
        cnt[f"{k}_{scope}"] = x[m].groupby(["panel_district_id", "year"]).size()
cnt = pd.DataFrame(cnt).fillna(0).reset_index().rename(columns={"panel_district_id": "district_id"})

d = load()
d = d.merge(cnt, on=["district_id", "year"], how="left")
for c in cnt.columns[2:]:
    d[c] = d[c].fillna(0)
    d.loc[d.year == 2025, c] = np.nan
print("\ncheck: 'polio_violent' equals the panel's polio_attacks (violent events only):",
      bool(np.allclose(d.polio_violent.fillna(-1), d.polio_attacks.fillna(-1))))
d["half"] = np.where(d.year <= 2021, "early", "late")

rows = []
print("\nTemporal ordering (Table S5 design) with each definition")
for expo, lab in [("polio_all", "polio-campaign, all event types (as Table S5)"),
                  ("polio_violent_or_riot", "polio-campaign, violent events or riots"),
                  ("polio_violent", "polio-campaign, violent events only"),
                  ("hw_all", "health-worker, all event types (approx. Table S5)"),
                  ("hw_violent_or_riot", "health-worker, violent events or riots"),
                  ("hw_violent", "health-worker, violent events only")]:
    tot = d.groupby("half")[expo].sum()
    for oh in ["late", "early"]:
        eh = "early" if oh == "late" else "late"
        xx = d[d.half == eh].groupby("district_id")[expo].mean()
        s_ = d[d.half == oh].copy()
        s_["x"] = s_.district_id.map(xx)
        s_ = s_.dropna(subset=["x", "ghsl_population"])
        s_["lpop"] = np.log(s_.ghsl_population)
        kp, ky = s_.groupby("province").wpv1_cases.sum(), s_.groupby("year").wpv1_cases.sum()
        s_ = s_[s_.province.isin(kp[kp > 0].index) & s_.year.isin(ky[ky > 0].index)].copy()
        s_["prov"] = s_.province.astype(str).astype("category")
        s_["yr"] = s_.year.astype(str).astype("category")
        m = pois("wpv1_cases ~ x + C(prov) + C(yr)", s_)
        b, se = m.params["x"], m.bse["x"]
        r = (np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se), m.pvalues["x"])
        print(line(f"{lab[:38]}: {eh}->{oh}", *r))
        rows.append(dict(exposure=lab, direction=f"{eh} exposure -> {oh} cases",
                         events_early=tot["early"], events_late=tot["late"],
                         irr=r[0], lo=r[1], hi=r[2], p=r[3]))
pd.DataFrame(rows).to_csv(f"{TAB}/review3_A2_temporal_order_by_definition.csv", index=False)
