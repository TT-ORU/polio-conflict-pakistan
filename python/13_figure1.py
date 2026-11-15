"""13  Figure 1: national WPV1 cases (bars) and ACLED violent events (line), from
outputs/tables/table1.csv, so each ACLED event is counted once. The earlier
figure summed the panel and so counted Karachi's city-wide events five times."""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import TAB, FIG

t = pd.read_csv(f"{TAB}/table1.csv", index_col=0)
fig, ax = plt.subplots(figsize=(7.3, 4))
red, navy = "#c0392b", "#2c3e50"
b = ax.bar(t.index, t.wpv1_cases, color=red, width=0.6)
for x, y in zip(t.index, t.wpv1_cases):
    ax.text(x, y + 2, f"{int(y)}", ha="center", va="bottom", fontsize=8, color="#7b241c")
ax.set_xlabel("Year")
ax.set_ylabel("Reported WPV1 cases", color=red)
ax.tick_params(axis="y", colors=red)
ax.set_xticks(t.index)
ax.set_ylim(0, 165)
ax2 = ax.twinx()
v = t.violent_events.dropna()
ax2.plot(v.index, v.values, color=navy, marker="o", lw=2)
ax2.set_ylabel("ACLED violent events (to 2024)", color=navy)
ax2.tick_params(axis="y", colors=navy)
ax2.set_ylim(0, 3000)
ax.spines["top"].set_visible(False); ax2.spines["top"].set_visible(False)
ax.set_title("National WPV1 cases and violent conflict, 2018–2025", fontsize=11)
fig.tight_layout()
fig.savefig(f"{FIG}/Figure1_national_trend.png", dpi=300)
fig.savefig(f"{FIG}/Figure1_national_trend.pdf")
print("wrote Figure 1")
