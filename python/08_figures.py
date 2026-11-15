"""08  Figure 3 (interaction: marginal effect and tertile grid) and Figure 4
(district map). Figures 1 and 2 are drawn by R/04_figures.R.

Figure 4 palette checked with a colour-vision validator: the three category
colours pass lightness, chroma, colour-blind separation and contrast."""
import numpy as np, pandas as pd, geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, NullFormatter
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from common import (load, estimation_sample, prep, pois, FIG, DATA, F_INT,
                    classify_temporal)

INK, NEITHER, VIOL, COV, BOTH, OUTSIDE = ("#1A1A1A", "#E8E6E1", "#2E6FAF",
                                          "#D97706", "#9B1B3C", "#F5F4F2")
d = load()
E = prep(estimation_sample(d))
m = pois(F_INT, E)
V = m.cov_params()
dist = E.groupby("district_id").agg(vm=("vm", "first"), zd=("zero_dose_rate_mics", "first"),
                                    cases=("wpv1_cases", "sum"), py=("ghsl_population", "sum"),
                                    name=("district_name", "first"))
dist["vt"] = pd.qcut(dist.vm, 3, labels=["low", "mid", "high"])
dist["zt"] = pd.qcut(dist.zd, 3, labels=["low", "mid", "high"])


def marg(zpct):
    z = zpct / 10
    b = m.params["vm10"] + z * m.params["vm10:zd10"]
    se = np.sqrt(V.loc["vm10", "vm10"] + z ** 2 * V.loc["vm10:zd10", "vm10:zd10"]
                 + 2 * z * V.loc["vm10", "vm10:zd10"])
    return np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se)


# ---- Figure 3 ---------------------------------------------------------------
zz = np.linspace(dist.zd.min() * 100, dist.zd.max() * 100, 200)
ir = np.array([marg(z) for z in zz])
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.2, 4.5),
                             gridspec_kw=dict(width_ratios=[1.35, 1]))
a1.fill_between(zz, ir[:, 1], ir[:, 2], color="#4477aa", alpha=.18, lw=0)
a1.plot(zz, ir[:, 0], color="#1f4e79", lw=2)
a1.axhline(1, color="grey", lw=.9, ls="--")
a1.set_yscale("log")
a1.set_ylim(0.42, 2.4)
a1.set_yticks([0.5, 0.75, 1, 1.5, 2])
a1.yaxis.set_major_formatter(ScalarFormatter())
a1.yaxis.set_minor_formatter(NullFormatter())
a1.plot(dist.zd * 100, np.full(len(dist), 0.455), "|", color="#333333", ms=7, mew=.9, alpha=.55)
a1.set_xlabel("District zero-dose rate (%)")
a1.set_ylabel("IRR per +10 violent events per year")
a1.set_title("A. Conflict association by immunisation coverage", loc="left", fontsize=10.5)
a1.text(.98, .10, "Band: 95% CI. Ticks: observed districts.", transform=a1.transAxes,
        ha="right", fontsize=8, color="#555555")
for sp in ["top", "right"]:
    a1.spines[sp].set_visible(False)
o = ["low", "mid", "high"]
M = np.zeros((3, 3)); C = np.zeros((3, 3), int); N = np.zeros((3, 3), int)
for i, vt in enumerate(o):
    for j, zt in enumerate(o):
        s = dist[(dist.vt == vt) & (dist.zt == zt)]
        M[i, j], C[i, j], N[i, j] = s.cases.sum() / s.py.sum() * 1e6, s.cases.sum(), len(s)
im = a2.imshow(M, cmap="YlOrRd", origin="lower", aspect="auto")
for i in range(3):
    for j in range(3):
        a2.text(j, i, f"{M[i,j]:.1f}\n{C[i,j]} cases\n{N[i,j]} distr.", ha="center",
                va="center", fontsize=8.5, color="white" if M[i, j] > M.max() * .6 else "#222222")
a2.set_xticks(range(3), o); a2.set_yticks(range(3), o)
a2.set_xlabel("Zero-dose rate tertile"); a2.set_ylabel("Violent events tertile")
a2.set_title("B. WPV1 cases per million person-years", loc="left", fontsize=10.5)
fig.colorbar(im, ax=a2, fraction=.046, pad=.04)
fig.tight_layout()
fig.savefig(f"{FIG}/Figure3_conflict_coverage.png", dpi=400)
fig.savefig(f"{FIG}/Figure3_conflict_coverage.pdf")
plt.close(fig)

# ---- Figure 4 ---------------------------------------------------------------
g = gpd.read_file(f"{DATA}/pakistan_districts_gadm_level3.gpkg")
dist["rate"] = dist.cases / dist.py * 1e6
dist["cls"] = np.where((dist.vt == "high") & (dist.zt == "high"), "both",
              np.where(dist.vt == "high", "viol", np.where(dist.zt == "high", "cov", "neither")))
flag = classify_temporal(d, dist.index)
flag = flag.index[flag.high_risk == 1]
_cls = pd.read_csv("../outputs/tables/tableS7_classification_by_district.csv")      # written by 05
_pan = pd.read_csv("../data/polio_panel_2018_2025.csv")
_later = _pan[_pan.year.between(2022, 2025) & _pan.district_id.isin(_cls.district_id)]
TVC = int(_later[_later.district_id.isin(_cls[_cls.high_risk == 1].district_id)].wpv1_cases.sum())
TVT = int(_later.wpv1_cases.sum())
g = g.merge(dist[["rate", "cls"]], left_on="district_id", right_index=True, how="left")
g["flagged"] = g.district_id.isin(flag)
prov = g.dissolve(by="province")
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.4, 6.2))
g[g.cls.isna()].plot(ax=a1, color=OUTSIDE, edgecolor="white", linewidth=.25)
for k, c in {"neither": NEITHER, "viol": VIOL, "cov": COV, "both": BOTH}.items():
    g[g.cls == k].plot(ax=a1, color=c, edgecolor="white", linewidth=.25)
g[g.flagged].boundary.plot(ax=a1, color=INK, linewidth=1.1)
prov.boundary.plot(ax=a1, color="#6B6B6B", linewidth=.6)
a1.set_title("A. District classification by chronic violence\nand routine-immunisation failure",
             loc="left", fontsize=11, color=INK)
cnt = g.cls.value_counts()
a1.legend(handles=[
    Patch(facecolor=BOTH, label=f"Top tertile of both ({cnt.get('both', 0)})"),
    Patch(facecolor=VIOL, label=f"Top tertile of violence only ({cnt.get('viol', 0)})"),
    Patch(facecolor=COV, label=f"Top tertile of zero-dose only ({cnt.get('cov', 0)})"),
    Patch(facecolor=NEITHER, label=f"Neither ({cnt.get('neither', 0)})"),
    Patch(facecolor=OUTSIDE, label=f"Not in estimation sample ({int(g.cls.isna().sum())})"),
    Line2D([0], [0], color=INK, lw=1.1, label=f"Flagged from 2018–2021 data only ({int(g.flagged.sum())})")],
    loc="upper center", bbox_to_anchor=(.5, -.01), ncol=2, fontsize=8, frameon=False,
    handlelength=1.5, columnspacing=1.6)
for txt, xy, xyt in [("Quetta block", (66.9, 30.1), (61.4, 26.4)),
                     ("Waziristan and\nsouthern KP", (70.1, 32.6), (74.6, 35.4))]:
    a1.annotate(txt, xy=xy, xytext=xyt, fontsize=8.2, color=INK, ha="center", va="center",
                arrowprops=dict(arrowstyle="-", color="#555555", lw=.7))
seq = LinearSegmentedColormap.from_list("seq", ["#FBF2F4", BOTH])
vmax = float(np.nanpercentile(g.rate, 98))
g[g.rate.isna()].plot(ax=a2, color=OUTSIDE, edgecolor="white", linewidth=.25)
g[g.rate.notna()].plot(ax=a2, column="rate", cmap=seq, edgecolor="white", linewidth=.25,
                       vmin=0, vmax=vmax)
prov.boundary.plot(ax=a2, color="#6B6B6B", linewidth=.6)
a2.set_title("B. Reported WPV1 incidence, 2019–2025\n(cases per million person-years)",
             loc="left", fontsize=11, color=INK)
cb = fig.colorbar(plt.cm.ScalarMappable(cmap=seq, norm=Normalize(0, vmax)), ax=a2,
                  fraction=.035, pad=.02)
cb.outline.set_visible(False)
cb.ax.tick_params(labelsize=8, length=0)
for ax in (a1, a2):
    ax.set_axis_off(); ax.set_aspect("equal")
fig.text(.012, .012, "Tertiles of district mean prior-year violent events and of the baseline "
         "zero-dose rate, among the 120 districts in the estimation sample. Dark outline in A marks "
         f"the {len(flag)} districts\nidentified from 2018–2021 information alone, which carried {TVC} of the "
         f"{TVT} cases reported in 2022–2025. Province boundaries in dark grey.",
         fontsize=7.4, color="#555555")
fig.tight_layout(rect=[0, .105, 1, 1])
fig.savefig(f"{FIG}/Figure4_map.png", dpi=400)
fig.savefig(f"{FIG}/Figure4_map.pdf")
print("wrote Figure3_conflict_coverage and Figure4_map (.png, .pdf)")
