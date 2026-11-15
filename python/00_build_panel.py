"""00  Build the analysis panel, 2018-2025, from the source files in data/.

This is a line-for-line port of R/01_consolidate_cases.R and R/02_build_panel.R,
extended to 2025. Those R scripts stop at 2024; the 2025 extension used in the
manuscript was built outside them, so this script makes that step explicit.

Two 2025 inputs are NOT derivable from the other source files:
  - 2025 district case counts: data/polio_cases_2025_district.csv
    (31 cases in 14 districts; extracted from the frozen panel. The author
    should record where these came from.)
  - 2025 population: the 2024 GHSL value carried forward unchanged, because
    data/district_population_ghsl.csv ends in 2024.
ACLED ends on 31 December 2024, so 2025 conflict is missing (not zero), and the
2025 regression year enters through its 2024 lag.

Writes data/polio_panel_2018_2025.csv and checks it cell by cell against the
frozen panel used for the manuscript.
"""
import sys
import numpy as np
import pandas as pd
from common import DATA

YEARS = range(2018, 2026)
KNOWN = {2018: 12, 2019: 147, 2020: 84, 2021: 1, 2022: 20, 2023: 6,
         2024: 74, 2025: 31}

frame = pd.read_csv(f"{DATA}/district_panel_frame.csv")
frame["district_id"] = frame.district_id.astype(int)
xw = pd.read_csv(f"{DATA}/polio_name_crosswalk.csv")
xw["polio_district"] = xw.polio_district.str.split().str.join(" ")

# ---- cases (R/01) ----------------------------------------------------------
# 2018-2021: weekly reports are cumulative, so take the year-end MAX, not a sum
raw = pd.read_csv(f"{DATA}/polio_raw_parsed.csv")
raw["district"] = raw.district.str.split().str.join(" ")
c1 = (raw[(raw.year >= 2018) & (raw.year <= 2021)]
      .groupby(["year", "district"])["count"].max().rename("cases").reset_index())
# 2022-2024: annual line list
ll = pd.read_csv(f"{DATA}/polio_cases_district_year.csv")
ll["district"] = ll.district.str.split().str.join(" ")
c2 = ll[(ll.year >= 2022) & (ll.year <= 2024)][["year", "district", "cases"]]
cases = pd.concat([c1, c2])
cases = cases[cases.cases > 0]
miss = set(cases.district) - set(xw.polio_district)
if miss:
    sys.exit(f"case districts missing from crosswalk: {sorted(miss)}")
cases = (cases.merge(xw, left_on="district", right_on="polio_district")
         .groupby(["panel_district_id", "year"]).cases.sum()
         .rename("wpv1_cases").reset_index()
         .rename(columns={"panel_district_id": "district_id"}))
# 2025: explicit input (see docstring)
c25 = pd.read_csv(f"{DATA}/polio_cases_2025_district.csv")
cases = pd.concat([cases, c25[["district_id", "year", "cases"]]
                   .rename(columns={"cases": "wpv1_cases"})])

panel = frame.merge(pd.DataFrame({"year": list(YEARS)}), how="cross")
panel = panel.merge(cases, on=["district_id", "year"], how="left")
panel["wpv1_cases"] = panel.wpv1_cases.fillna(0).astype(int)

tot = panel.groupby("year").wpv1_cases.sum()
bad = {y: (int(tot[y]), k) for y, k in KNOWN.items() if tot[y] != k}
if bad:
    sys.exit(f"national checksum failed (panel, known): {bad}")
print("national checksum OK:", tot.to_dict())

# ---- exposures, population, coverage (R/02) --------------------------------
ax = pd.read_csv(f"{DATA}/acled_to_panel_district_crosswalk.csv")
ax["acled_district"] = ax.acled_district.str.split().str.join(" ")
gen = pd.read_csv(f"{DATA}/acled_pakistan_district_year.csv")
gen["district"] = gen.district.str.split().str.join(" ")
gen["conflict_viol"] = gen.n_battle + gen.n_vac + gen.n_explo
gen = (gen.merge(ax, left_on="district", right_on="acled_district")
       .groupby(["panel_district_id", "year"])
       .agg(conflict_events=("n_events", "sum"),
            conflict_viol=("conflict_viol", "sum"),
            conflict_fatal=("fatalities", "sum")).reset_index()
       .rename(columns={"panel_district_id": "district_id"}))
# attacks: violent events only (built by 00a_build_attack_counts.py). The frozen
# panel used keyword matches over all event types, mostly protests.
pol = pd.read_csv(f"{DATA}/acled_attacks_district_year.csv")
pol["district"] = pol.district.str.split().str.join(" ")
pol = (pol.merge(ax, left_on="district", right_on="acled_district")
       .groupby(["panel_district_id", "year"])
       .agg(polio_attacks=("n_polio", "sum"),
            hworker_attacks=("n_hworker", "sum"),
            polio_attack_deaths=("polio_deaths", "sum")).reset_index()
       .rename(columns={"panel_district_id": "district_id"}))
pop = pd.read_csv(f"{DATA}/district_population_ghsl.csv")
pop25 = pop[pop.year == 2024].assign(year=2025)          # carried forward
pop = pd.concat([pop, pop25])
cov = pd.read_csv(f"{DATA}/district_zero_dose_mics.csv")   # corrected rate, children 12-23 months (review4)

panel = (panel.merge(gen, on=["district_id", "year"], how="left")
         .merge(pol, on=["district_id", "year"], how="left")
         .merge(pop, on=["district_id", "year"], how="left")
         .merge(cov, on="district_id", how="left"))
conf = ["conflict_events", "conflict_viol", "conflict_fatal",
        "polio_attacks", "hworker_attacks", "polio_attack_deaths"]
for c in conf:
    panel[c] = panel[c].fillna(0)
    panel.loc[panel.year == 2025, c] = np.nan            # ACLED ends 2024
panel = panel.sort_values(["district_id", "year"]).reset_index(drop=True)
g = panel.groupby("district_id")
panel["conflict_viol_lag1"] = g.conflict_viol.shift(1)
panel["conflict_fatal_lag1"] = g.conflict_fatal.shift(1)
panel["hworker_attacks_lag1"] = g.hworker_attacks.shift(1)
panel["hworker_attacks_lag2"] = g.hworker_attacks.shift(2)
panel["polio_attacks_lag1"] = g.polio_attacks.shift(1)
panel["lpop"] = np.log(panel.ghsl_population)

frozen = pd.read_csv(f"{DATA}/polio_panel_2018_2025_FROZEN.csv")
panel = panel[frozen.columns]
panel.to_csv(f"{DATA}/polio_panel_2018_2025.csv", index=False)

# ---- does the rebuild reproduce the frozen panel? --------------------------
a = panel.sort_values(["district_id", "year"]).reset_index(drop=True)
b = frozen.sort_values(["district_id", "year"]).reset_index(drop=True)
diffs = []
for c in b.columns:
    if pd.api.types.is_numeric_dtype(b[c]):
        same = np.isclose(a[c], b[c], rtol=1e-9, atol=1e-9, equal_nan=True)
    else:
        same = (a[c].astype(str) == b[c].astype(str)).values
    if not same.all():
        diffs.append((c, int((~same).sum())))
print(f"rebuilt panel: {a.shape[0]} rows x {a.shape[1]} columns")

# hworker_attacks_lag2 differed from the frozen panel in 2025 (a three-year lag there); it
# is now one of the redefined attack columns below, and no model uses it.
# The attack columns were redefined in the third revision (violent events
# only), so they are expected to differ from the frozen panel.
ATTACKS = {"polio_attacks", "hworker_attacks", "polio_attack_deaths",
           "hworker_attacks_lag1", "hworker_attacks_lag2", "polio_attacks_lag1"}
# The zero-dose rate was corrected in the fourth revision (review4/PRESPECIFICATION_zero_dose_correction.md):
# never-vaccinated children are now counted, and the age group is 12-23 months. It is expected to differ.
ATTACKS = ATTACKS | {"zero_dose_rate_mics"}
redefined = [d for d in diffs if d[0] in ATTACKS]
if redefined:
    print("  redefined attack columns (violent events only):", redefined)
other = [d for d in diffs if d[0] not in ATTACKS]
if other:
    print("UNEXPECTED differences from the frozen panel:", other)
    sys.exit(1)
print("rebuilt panel reproduces the frozen manuscript panel in every other "
      "cell of every column")
