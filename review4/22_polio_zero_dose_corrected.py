# Record copy. Run from the zero-dose publication package (code/), which holds the MICS-derived inputs.
"""22  Corrected district zero-dose rates for the companion polio paper (polio review4 specification, sha256 cc565d70...).
Corrected definition (never-vaccinated reports counted as zero-dose). Primary 12-23 months; sensitivity 12-35 and
0-35 months. MICS districts mapped to polio panel districts via master_district_ids. District aggregates only."""
import numpy as np, pandas as pd
from paths import INPUTS as U, RESTRICTED
ch = pd.read_pickle("out/children_exposure.pkl")
dist = pd.read_csv("out/aim1_district.csv")[["dkey", "master_district_ids"]]
link = dist.assign(district_id=dist.master_district_ids.astype(str).str.split(";")).explode("district_id")
link["district_id"] = link.district_id.astype(int)
C = ch[ch.zero_dose.notna()][["dkey", "cage", "zero_dose", "chweight"]].merge(link[["dkey", "district_id"]], on="dkey")
out = None
for lab, lo, hi in [("1223", 12, 23), ("1235", 12, 35), ("0035", 0, 35)]:
    g = C[C.cage.between(lo, hi)].groupby("district_id")
    t = pd.DataFrame({f"zdw_{lab}": g.apply(lambda x: np.average(x.zero_dose, weights=x.chweight), include_groups=False),
                      f"zdu_{lab}": g.zero_dose.mean(), f"n_{lab}": g.size()})
    out = t if out is None else out.join(t, how="outer")
old = pd.read_csv(U + "polio/zero_dose_check_review3.csv").set_index("district_id")
out = out.join(old[["zdw"]].rename(columns={"zdw": "zdw_previous"}), how="outer")
print("panel districts:", len(out), " with corrected 12-23:", out.zdw_1223.notna().sum(), " previous:", out.zdw_previous.notna().sum())
print(out[["zdw_1223", "zdw_1235", "zdw_0035", "zdw_previous"]].describe().round(3).to_string())
print("correlations:\n", out[["zdw_1223", "zdw_1235", "zdw_0035", "zdw_previous"]].corr().round(3).to_string())
out.reset_index().to_csv("out/polio_zero_dose_check_review4.csv", index=False)
out[["zdw_1223"]].rename(columns={"zdw_1223": "zero_dose_rate_mics"}).dropna().round(4).reset_index().to_csv(
    "out/polio_district_zero_dose_mics_corrected.csv", index=False)
