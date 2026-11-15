"""07  Spatial sensitivity analysis (Table S9): neighbouring-district violence
and cases, and Moran's I on the residuals of the primary model, by year.
Queen contiguity over all 141 districts from data/pakistan_districts_gadm_level3.gpkg."""
import numpy as np, pandas as pd, geopandas as gpd
from libpysal.weights import Queen
from esda.moran import Moran
from common import (load, estimation_sample, prep, pois, line, Keys, TAB, DATA,
                    F_INT, SEED)

K = Keys("07_spatial")
g = gpd.read_file(f"{DATA}/pakistan_districts_gadm_level3.gpkg").sort_values(
    "district_id").reset_index(drop=True)
d = load()
assert set(g.district_id) == set(d.district_id), "boundary IDs do not match panel"
W = Queen.from_dataframe(g, use_index=False)
nb = {int(g.district_id[i]): [int(g.district_id[j]) for j in W.neighbors[i]]
      for i in range(len(g))}
n = pd.Series({k: len(v) for k, v in nb.items()})
print(f"Queen contiguity: median {n.median():.0f} neighbours "
      f"(range {n.min()}-{n.max()}), {(n == 0).sum()} islands")
K.add("nb_median", n.median())

viol = d.pivot_table(index="district_id", columns="year", values="conflict_viol_lag1")
cases = d.pivot_table(index="district_id", columns="year", values="wpv1_cases")


def neigh(tbl, i, yr):
    if yr not in tbl.columns:
        return np.nan
    v = [tbl.loc[j, yr] for j in nb[i] if j in tbl.index and pd.notna(tbl.loc[j, yr])]
    return np.mean(v) if v else np.nan


d["nb_viol_lag1"] = [neigh(viol, r.district_id, r.year) for r in d.itertuples()]
d["nb_cases_lag1"] = [neigh(cases, r.district_id, r.year - 1) for r in d.itertuples()]
E = prep(estimation_sample(d))
E["nm"] = E.groupby("district_id").nb_viol_lag1.transform("mean")
E["nm10"], E["nd10"] = E.nm / 10, (E.nb_viol_lag1 - E.nm) / 10
E["nbc"] = np.log1p(E.nb_cases_lag1)

rows = []
for key, lab, f, terms in [
        ("s0", "no neighbour terms", F_INT, ["vm10:zd10"]),
        ("s1", "with neighbour violence",
         "wpv1_cases ~ vm10*zd10 + vd10 + nm10 + nd10 + C(prov) + C(yr)",
         ["vm10:zd10", "nm10", "nd10"]),
        ("s2", "with neighbour violence and cases",
         "wpv1_cases ~ vm10*zd10 + vd10 + nm10 + nd10 + nbc + C(prov) + C(yr)",
         ["vm10:zd10", "nbc"])]:
    m = pois(f, E)
    print(f" {lab}:")
    for t in terms:
        r = K.ci(f"{key}_{t.replace(':', 'x')}", m, t)
        print(line(t, *r))
        rows.append(dict(model=lab, term=t, irr=r[0], lo=r[1], hi=r[2], p=r[3]))
pd.DataFrame(rows).to_csv(f"{TAB}/tableS9_spatial.csv", index=False)

m0 = pois(F_INT, E)
gi = g.set_index("district_id")
mor = []
np.random.seed(SEED)
print(" Moran's I, Pearson residuals of the primary model:")
for yr, grp in E.groupby("year"):
    mu = m0.predict(grp, offset=grp.lpop)
    r = (grp.wpv1_cases - mu) / np.sqrt(mu.clip(lower=1e-9))
    w = Queen.from_dataframe(gi.loc[grp.district_id].reset_index(), use_index=False)
    w.transform = "r"
    mi = Moran(r.values, w, permutations=999)
    mor.append(dict(year=int(yr), I=mi.I, p=mi.p_sim))
    K.add(f"moran_{int(yr)}", mi.I)
    print(f"  {int(yr)}  I = {mi.I:+.3f}  p = {mi.p_sim:.3f}")
mor = pd.DataFrame(mor)
K.add("moran_mean", mor.I.mean())
K.add("moran_n_sig", (mor.p < 0.05).sum())
mor.to_csv(f"{TAB}/tableS9_moran_by_year.csv", index=False)
print(f"  mean I = {mor.I.mean():+.3f}; significant in {(mor.p < .05).sum()} of 7 years")
K.save()
