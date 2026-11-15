"""12  Table 1: national WPV1 cases and ACLED indicators by year.

Each ACLED event is counted once. (Summing the panel would count Karachi's
city-wide events once for each of its five districts.) Only ACLED districts
that map to the panel geography are counted, as in the models.
"""
import pandas as pd
from common import load, DATA, TAB

d = load()
ax = pd.read_csv(f"{DATA}/acled_to_panel_district_crosswalk.csv")
ax["acled_district"] = ax.acled_district.str.split().str.join(" ")
mapped = set(ax.acled_district)

gen = pd.read_csv(f"{DATA}/acled_pakistan_district_year.csv")
gen["district"] = gen.district.str.split().str.join(" ")
gen = gen[gen.district.isin(mapped)]
gen["viol"] = gen.n_battle + gen.n_vac + gen.n_explo
att = pd.read_csv(f"{DATA}/acled_attacks_district_year.csv")
att = att[att.district.isin(mapped)]

t = pd.DataFrame({"wpv1_cases": d.groupby("year").wpv1_cases.sum()})
t["violent_events"] = gen.groupby("year").viol.sum()
t["health_worker_attacks"] = att.groupby("year").n_hworker.sum()
t["polio_campaign_attacks"] = att.groupby("year").n_polio.sum()
t = t.loc[2018:2025]
t.loc[2025, ["violent_events", "health_worker_attacks", "polio_campaign_attacks"]] = None
t.to_csv(f"{TAB}/table1.csv")
print(t.to_string())
from common import Keys
K = Keys("12_table1")
for y, r in t.iterrows():
    for c, k in [("violent_events", "viol"), ("health_worker_attacks", "hw"), ("polio_campaign_attacks", "pol")]:
        if pd.notna(r[c]):
            K.add(f"t1_{k}_{y}", r[c])
K.save()
