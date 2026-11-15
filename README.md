# Armed conflict, routine immunisation gaps and reported wild poliovirus incidence in Pakistan

Code and district-level data for the article *Armed conflict, routine immunisation gaps and reported wild
poliovirus incidence in Pakistan: a district-level panel analysis, 2018-2025* by Faisal Munir, Tor Arnison,
Syed Taha Jamil Abidi, Jahangir Khan and Tahir Taj (reference to follow).

The pipeline builds a panel of 141 districts for 2018-2025, fits Poisson within-between models of reported
wild poliovirus type 1 (WPV1) cases on lagged ACLED violent events and MICS zero-dose prevalence, and
produces every table and figure. `check_against_frozen.py` compares a fresh run with the frozen results, and
`python/99_verify.py` checks every number reported in the article against the run.

## Contents

| Path | What it holds |
|---|---|
| `python/` | Analysis scripts, run in order by `run_all.sh`; `manuscript_numbers.csv` lists every number reported in the article |
| `python/get_boundaries.py` | Downloads the GADM district boundaries, which may not be redistributed (see below) |
| `data/` | District-year inputs (see "Data sources and terms") |
| `frozen_results/` | Tables, keys, figures and logs of the frozen run of 21 October 2026 |
| `review3/`, `review4/` | Analyses specified in writing before they were run, with timestamps and SHA-256; the script that computed district zero-dose rates from MICS microdata |
| `environment/` | Python package versions; SHA-256 of the raw ACLED extract held by the authors |
| `docs/` | Notes from the frozen package: revision history, reproduction check, original manifest |
| `run_all.sh`, `check_against_frozen.py` | Run everything; compare the run with the frozen results |

## How to run

Python 3.13 with the package versions in `requirements.txt`.

```
pip install -r requirements.txt
python3 python/get_boundaries.py    # once: downloads GADM 4.1 boundaries for Pakistan
bash run_all.sh                     # about 10 minutes; FAST=1 bash run_all.sh takes about 1 minute
python3 check_against_frozen.py     # every table, key and figure against frozen_results/
```

Results are written to `outputs/`. The randomisation inference and bootstraps use fixed seeds, so a full run
reproduces the frozen results exactly. With `FAST=1` the randomisation inference uses fewer permutations, so
its results differ by design.

Two steps need the raw ACLED event file, which may not be redistributed: `python/00a_build_attack_counts.py`
(its output, `data/acled_attacks_district_year.csv`, is included) and `python/11_review3_acled_attacks.py`.
Researchers with their own ACLED access can rerun them; the SHA-256 of the authors' extract is in
`environment/restricted_inputs.sha256`.

## Data sources and terms

All files in `data/` are district or district-year aggregates. No individual-level records are included.

| Files | Source | Terms |
|---|---|---|
| `polio_raw_parsed.csv`, `polio_cases_district_year.csv`, `polio_cases_2025_district.csv`, `polio_name_crosswalk.csv` | Reported WPV1 cases by district and year, from WHO Eastern Mediterranean Regional Office polio bulletins and district case lists | Public sources; cite WHO |
| `acled_pakistan_district_year.csv`, `acled_polio_district_year.csv`, `acled_attacks_district_year.csv`, `acled_to_panel_district_crosswalk.csv` | Event counts by district and year, aggregated by the authors from the Armed Conflict Location & Event Data project (ACLED) | ACLED Terms & Conditions; raw events are not included. Cite ACLED as its attribution policy requires (www.acleddata.com) |
| `district_zero_dose_mics.csv`, `district_zero_dose_mics_v03_original.csv`, `zero_dose_check_review3.csv`, `zero_dose_check_review4.csv` | District zero-dose rates computed by the authors from UNICEF Multiple Indicator Cluster Surveys (MICS) | MICS microdata are available from UNICEF after registration (https://mics.unicef.org/surveys) and are not included |
| `district_population_ghsl.csv` | District population from GHS-POP R2023A (Schiavina et al. 2023, European Commission Joint Research Centre) | CC BY 4.0 |
| `district_covariates_review3.csv` | Travel time to the nearest city (2015), OpenStreetMap health facilities and VIIRS night-time lights, by district | Cite each source. Health-facility counts are derived from OpenStreetMap (© OpenStreetMap contributors) and are available under the Open Database License |
| `district_panel_frame.csv`, `polio_panel_2018_2025.csv`, `polio_panel_2018_2025_FROZEN.csv` | Analysis panel built by `python/00_build_panel.py` from the files above | As for the files above |

## Not included

- Raw ACLED event data, which ACLED's terms do not allow to be redistributed.
- GADM district boundaries (https://gadm.org/license.html); `python/get_boundaries.py` rebuilds the file.
- MICS microdata (UNICEF terms).
- Manuscript drafts and correspondence with the journal.

## Changes from the frozen package (version 0.6, 21 October 2026)

- The items listed under "Not included" were removed, together with a script that edited the manuscript
  (`python/patch5_corrections.py`), and review reports and scripts for manuscript revision.
- `python/14_figure2.py` was added (Figure 2, 14 November 2026: districts with no cases or no violent events are
  plotted at zero) and is run by `run_all.sh`. It has no frozen counterpart, so `check_against_frozen.py` does
  not check it.
- `run_all.sh` stops with a message if the boundary file is missing.
- The original manifest is `docs/FROZEN_MANIFEST_v0.6.sha256` (SHA-256 of the manifest:
  ef2fc425908125c6824eb7d5e4f9d358104a5d462b687c61159d4cb156c217e0). Files carried over unchanged can be
  checked against it; files that were removed are listed in it but are not in this repository.

## Licence

Code: MIT (see `LICENSE`). Data: the terms of each source, as listed above.

## Contact

Tahir Taj, Department of Clinical Epidemiology and Biostatistics, Faculty of Medicine and Health, Örebro
University, Sweden. tahir.taj@oru.se
