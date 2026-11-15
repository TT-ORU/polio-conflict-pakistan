"""Shared paths, sample definition and model helpers for every analysis script.

Change something here and it changes everywhere, which is the point: the
sample, the Mundlak split and the fixed-effects handling are defined once.
"""
import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
TAB = os.path.join(ROOT, "outputs", "tables")
FIG = os.path.join(ROOT, "outputs", "figures")
KEYS = os.path.join(ROOT, "outputs", "keys")
for p in (TAB, FIG, KEYS):
    os.makedirs(p, exist_ok=True)

SEED = 20261016
FAST = os.environ.get("FAST") == "1"     # FAST=1 shortens the resampling runs

F_MAIN = "wpv1_cases ~ vm10 + vd10 + zd10 + C(prov) + C(yr)"
F_INT = "wpv1_cases ~ vm10*zd10 + vd10 + C(prov) + C(yr)"


def load():
    return pd.read_csv(os.path.join(DATA, "polio_panel_2018_2025.csv"))


def estimation_sample(d):
    """2019-2025, lagged violence, population and MICS zero-dose all present.
    840 district-years, 120 districts, 361 cases."""
    e = d[(d.year >= 2019) & d.conflict_viol_lag1.notna()
          & d.ghsl_population.notna() & d.zero_dose_rate_mics.notna()].copy()
    assert len(e) == 840 and e.district_id.nunique() == 120, "sample changed"
    return e.sort_values(["district_id", "year"])


def prep(s, var="conflict_viol_lag1", drop_empty_prov=True):
    """Mundlak split of lagged violence within the sample given. vm10/vd10 are
    per 10 events/yr, vm1/vd1 per event; zd10 is per 10 percentage points.
    Province and year levels with no cases are dropped: their fixed effect is
    minus infinity and the model would otherwise fail to converge. Pass
    drop_empty_prov=False for a model with no province effects."""
    s = s.copy()
    s["vm"] = s.groupby("district_id")[var].transform("mean")
    s["vm1"] = s.vm
    s["vd1"] = s[var] - s.vm
    s["vm10"] = s.vm1 / 10
    s["vd10"] = s.vd1 / 10
    if "zero_dose_rate_mics" in s:
        s["zd10"] = s.zero_dose_rate_mics * 10
    s["lpop"] = np.log(s.ghsl_population)
    kp = s.groupby("province").wpv1_cases.sum()
    ky = s.groupby("year").wpv1_cases.sum()
    keep = s.year.isin(ky[ky > 0].index)
    if drop_empty_prov:
        keep &= s.province.isin(kp[kp > 0].index)
    s = s[keep].copy()
    s["prov"] = s.province.astype(str).astype("category")
    s["yr"] = s.year.astype(str).astype("category")
    return s


def pois(formula, data, cluster=True):
    kw = (dict(cov_type="cluster", cov_kwds={"groups": data.district_id})
          if cluster else {})
    return smf.glm(formula, data=data, family=sm.families.Poisson(),
                   offset=data.lpop).fit(**kw)


def ci(m, term):
    b, se = m.params[term], m.bse[term]
    return np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se), m.pvalues[term]


def line(label, irr, lo, hi, p):
    return f"  {label:<56} {irr:7.3f} ({lo:.3f} to {hi:.3f})  p = {p:.4f}"


class Keys:
    """Collects named results; 99_verify.py compares them with the numbers
    printed in the manuscript."""

    def __init__(self, name):
        self.name, self.rows = name, []

    def add(self, key, value):
        self.rows.append((key, float(value)))

    def ci(self, key, m, term):
        irr, lo, hi, p = ci(m, term)
        for k, v in [("", irr), ("_lo", lo), ("_hi", hi), ("_p", p)]:
            self.add(key + k, v)
        return irr, lo, hi, p

    def save(self):
        pd.DataFrame(self.rows, columns=["key", "value"]).to_csv(
            os.path.join(KEYS, f"{self.name}.csv"), index=False)


def classify_temporal(d, ids):
    """High-risk flag from 2018-2021 violence and baseline zero-dose ONLY.
    Tertile cut-points across the given analysis districts; no case data."""
    tr = d[(d.year <= 2021) & d.conflict_viol.notna()]
    b = d.groupby("district_id").agg(zd=("zero_dose_rate_mics", "first"),
                                     name=("district_name", "first"))
    b["vm_train"] = tr.groupby("district_id").conflict_viol.mean()
    b = b.dropna(subset=["zd", "vm_train"])
    b = b[b.index.isin(ids)].copy()
    b["vt"] = pd.qcut(b.vm_train, 3, labels=["low", "mid", "high"])
    b["zt"] = pd.qcut(b.zd, 3, labels=["low", "mid", "high"])
    b["high_risk"] = ((b.vt == "high") & (b.zt == "high")).astype(int)
    return b
