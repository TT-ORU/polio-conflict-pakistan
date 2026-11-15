"""Rebuild data/pakistan_districts_gadm_level3.gpkg from GADM 4.1.

GADM boundaries may be used for academic work but may not be redistributed
(https://gadm.org/license.html), so the file is not in this repository.
This script downloads GADM 4.1 for Pakistan, keeps the 141 level-3 districts and
adds the panel's district_id by matching district names to
data/district_panel_frame.csv. It stops if any district does not match.

After `bash run_all.sh`, `python3 check_against_frozen.py` confirms that the
spatial tables and the Figure 4 map are identical to the frozen results.
"""
import os
import sys
import urllib.request
import pandas as pd
import geopandas as gpd

URL = "https://geodata.ucdavis.edu/gadm/gadm4.1/gpkg/gadm41_PAK.gpkg"
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
RAW = os.path.join(DATA, "gadm41_PAK.gpkg")                    # GADM download (not to be shared)
OUT = os.path.join(DATA, "pakistan_districts_gadm_level3.gpkg")

# download once
if not os.path.exists(RAW):
    print("downloading", URL)
    urllib.request.urlretrieve(URL, RAW)

# GADM level-3 units (districts) and the panel's list of 141 districts
g = gpd.read_file(RAW, layer="ADM_ADM_3")
frame = pd.read_csv(os.path.join(DATA, "district_panel_frame.csv"))

# match on district name; each panel district must match exactly one GADM unit
g["district_name"] = g.NAME_3.str.strip()
m = frame.merge(g[["district_name", "NAME_2", "geometry"]], on="district_name", how="left")
missing = m[m.geometry.isna()].district_name.tolist()
if len(g) != 141 or len(m) != 141 or missing:
    print("GADM units:", len(g), "| matched rows:", len(m), "| panel districts without a match:", missing)
    sys.exit("The boundaries do not match the panel. Please contact the corresponding author.")

# same columns and order as the file used in the analysis
out = gpd.GeoDataFrame({"district_id": m.district_id, "province": m.province,
                        "division": m.NAME_2, "district": m.district_name},
                       geometry=m.geometry.values, crs=g.crs)
out = out.sort_values("district_id").reset_index(drop=True)
out.to_file(OUT, layer="pakistan_districts_gadm_level3", driver="GPKG")
print("wrote", OUT, "with", len(out), "districts")
