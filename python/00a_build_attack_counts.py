"""00a  District-year counts of attacks on polio activity and on health workers.

Run once, before 00_build_panel.py, with the ACLED event file (not distributed
with this package under ACLED's terms of use):
    python 00a_build_attack_counts.py /path/to/acled_pakistan_events_notes.csv

An event counts as an attack only if it is a violent event (battles, violence
against civilians, explosions/remote violence). The earlier version of these
counts used every event type, so most "attacks" were protests; see
review3/Polio_review3_results.docx, check A2.

  polio attack:          notes mention polio, vaccinat* or immuni*
  health-worker attack:  an actor is "Health workers", or the notes mention
                         health worker(s), lady health worker, vaccinat*,
                         doctor, nurse or paramedic

Writes data/acled_attacks_district_year.csv on the ACLED district names used by
data/acled_to_panel_district_crosswalk.csv.
"""
import sys
import pandas as pd
from common import DATA

VIOLENT = ["Battles", "Violence against civilians", "Explosions/Remote violence"]
ev = pd.read_csv(sys.argv[1], low_memory=False)
notes = ev.notes.fillna("").str.lower()
actors = (ev.actor1.fillna("") + "|" + ev.assoc_actor_1.fillna("") + "|"
          + ev.actor2.fillna("") + "|" + ev.assoc_actor_2.fillna("")).str.lower()
violent = ev.event_type.isin(VIOLENT)
ev["polio"] = violent & notes.str.contains("polio|vaccinat|immuni")
ev["hworker"] = violent & (actors.str.contains("health workers")
                           | notes.str.contains(r"health worker|lady health|vaccinat|doctor|nurse|paramedic"))
ev["polio_deaths"] = ev.fatalities.where(ev.polio, 0)
ev["district"] = ev.admin2.fillna("").str.upper().str.split().str.join(" ")
out = (ev[ev.polio | ev.hworker]
       .groupby(["district", "year"])
       .agg(n_polio=("polio", "sum"), n_hworker=("hworker", "sum"),
            polio_deaths=("polio_deaths", "sum")).reset_index())
out.to_csv(f"{DATA}/acled_attacks_district_year.csv", index=False)
print(f"wrote {len(out)} district-years; polio attacks {out.n_polio.sum()}, "
      f"health-worker attacks {out.n_hworker.sum()} (all years)")
print(out.groupby("year")[["n_polio", "n_hworker"]].sum().loc[2017:].to_string())
