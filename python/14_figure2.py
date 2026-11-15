"""14  Figure 2: district co-location of polio burden and conflict.

Added 14 November 2026 (not part of the frozen v0.6 run). Writes Figure 2 to
outputs/figures and the plotted value for every district to outputs/tables.

x: mean annual ACLED violent events per district, 2018-2024 (ACLED ends in 2024).
y: total reported WPV1 cases per district, 2018-2025.
Both axes are linear from 0 to 1 and logarithmic above 1, so districts with no
cases or no violent events are drawn at 0, not dropped or shifted.
Karachi's five units share the city-wide ACLED count, so they share one x value.
"""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import SymmetricalLogLocator, NullLocator

from common import load, FIG, TAB

d = load()

# violence is complete for 2018-2024 and missing for 2025; cases are complete
assert d[d.year <= 2024].conflict_viol.notna().all()
assert d[d.year == 2025].conflict_viol.isna().all()
assert d.wpv1_cases.sum() == 375

# one row per district (the mean skips the missing 2025 value, so it covers 2018-2024)
dist = d.groupby(["district_id", "district_name", "province"], as_index=False).agg(
    violent_events_mean=("conflict_viol", "mean"), wpv1_cases=("wpv1_cases", "sum"))
assert len(dist) == 141

# one colour per province, as in the earlier version of the figure; the three small grey
# groups also get their own marker so they can be told apart
names = {"Azad Kashmir": "Azad Kashmir", "Balochistan": "Balochistan",
         "Federally Administered Tribal Ar": "Former FATA", "Gilgit-Baltistan": "Gilgit-Baltistan",
         "Islamabad": "Islamabad", "Khyber-Pakhtunkhwa": "Khyber Pakhtunkhwa",
         "Punjab": "Punjab", "Sindh": "Sindh"}
dist["group"] = dist.province.map(names)
assert dist.group.notna().all()
dist.sort_values(["group", "district_name"]).to_csv(
    os.path.join(TAB, "figure2_district_values.csv"), index=False)

style = [("Azad Kashmir", "#a3b1ae", "s"), ("Balochistan", "#9b59b6", "o"),
         ("Former FATA", "#f0973c", "o"), ("Gilgit-Baltistan", "#c9ccce", "^"),
         ("Islamabad", "#7f8c8d", "D"), ("Khyber Pakhtunkhwa", "#cc5a50", "o"),
         ("Punjab", "#4cbb7c", "o"), ("Sindh", "#4a90c8", "o")]

fig, ax = plt.subplots(figsize=(7.3, 4.65))
for name, colour, marker in style:
    g = dist[dist.group == name]
    ax.scatter(g.violent_events_mean, g.wpv1_cases, s=40, c=colour, marker=marker,
               edgecolors="white", linewidths=0.6, alpha=0.9,
               label=f"{name} ({len(g)})", zorder=3)

# linear between 0 and 1, logarithmic above 1, on both axes
ax.set_xscale("symlog", linthresh=1, linscale=0.6)
ax.set_yscale("symlog", linthresh=1, linscale=0.35)   # counts are whole numbers, so 0-1 can be short
ax.set_xticks([0, 1, 10, 100], labels=["0", "1", "10", "100"])
ax.set_yticks([0, 1, 2, 5, 10, 20, 50], labels=["0", "1", "2", "5", "10", "20", "50"])
ax.xaxis.set_minor_locator(SymmetricalLogLocator(base=10, linthresh=1, subs=range(2, 10)))
ax.yaxis.set_minor_locator(NullLocator())
ax.set_xlim(-0.12, dist.violent_events_mean.max() * 1.4)
ax.set_ylim(-0.12, dist.wpv1_cases.max() * 1.4)

# the districts with no violence and no cases overlap at the origin, so say how many
n00 = int(((dist.wpv1_cases == 0) & (dist.violent_events_mean == 0)).sum())
ax.annotate(f"{n00} districts with no violent events and no cases",
            xy=(0, 0), xytext=(0.06, 0.45), fontsize=7.5, color="#52514e",
            arrowprops=dict(arrowstyle="-", color="#8a8984", lw=0.6))

ax.grid(True, which="major", color="#e6e5e1", lw=0.7, zorder=0)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.set_xlabel("Mean annual violent events per district, 2018–2024")
ax.set_ylabel("Total reported WPV1 cases, 2018–2025")
ax.set_title("District co-location of polio burden and conflict", fontsize=11)
ax.legend(loc="upper left", ncol=2, fontsize=8, frameon=False, title="Province (districts)",
          title_fontsize=8, alignment="left")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "Figure2_colocation.png"), dpi=300)
fig.savefig(os.path.join(FIG, "Figure2_colocation.pdf"))

# numbers used in the caption
print("districts", len(dist), "| no cases", int((dist.wpv1_cases == 0).sum()),
      "| no violent events", int((dist.violent_events_mean == 0).sum()),
      "| neither", int(((dist.wpv1_cases == 0) & (dist.violent_events_mean == 0)).sum()))
print(dist.groupby("group").agg(districts=("district_id", "size"), cases=("wpv1_cases", "sum")))
print("max x", dist.violent_events_mean.max(), "| max y", dist.wpv1_cases.max())
