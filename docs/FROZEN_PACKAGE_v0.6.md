# Publication package: Conflict marks the reservoir, not the moment

Frozen on 21 October 2026. Manuscript version 0.6 (corrected zero-dose rate), prepared for Conflict and Health.
This file describes the frozen package; `README.md` documents the pipeline and each revision in detail.

## Contents

| Folder | What it holds |
|---|---|
| `manuscript/` | Manuscript v0.6 clean (docx, pdf) and tracked; `history/` holds tracked versions 0.3, 0.4 and 0.5 |
| `submission/` | Cover letter v0.5 (clean, tracked), the zero-dose correction note (old and new results), the list of every v0.4 edit |
| `python/` | The analysis pipeline (`run_all.sh` runs it); `manuscript_numbers.csv` lists all 367 numbers printed in the manuscript |
| `data/` | All inputs except the raw ACLED event file: case reconstruction, ACLED district-year aggregates, crosswalks, GHSL population, district boundaries, the corrected district zero-dose rate (from the zero-dose package, step 22) and the previous rate |
| `outputs/` | Tables, keys, figures and logs of the frozen run (overwritten by a rerun) |
| `frozen_results/` | A read-only copy of `outputs/` for comparison after a rerun |
| `review3/`, `review4/` | Pre-specifications with timestamps and SHA-256, review checks, the manuscript revision scripts and `verify_tables_v05.py` |
| `environment/` | Python package versions; SHA-256 of the raw ACLED file |
| `MANIFEST.sha256` | SHA-256 of every file |

## Reproducing

```
bash run_all.sh                    # about 10 minutes; 20,000 randomisation permutations
python3 check_against_frozen.py    # every table, key and figure against frozen_results/
python3 review4/verify_tables_v05.py manuscript/Polio_conflict_manuscript_v0.6_CLEAN.docx
```

`run_all.sh` ends with `python/99_verify.py`, which checks all 367 manuscript numbers against the run. The
randomisation inference and bootstraps use fixed seeds, so a rerun reproduces the frozen outputs exactly.

Two steps need the raw ACLED event file, which is not redistributed (ACLED terms): `python/00a_build_attack_counts.py`
(its output, `data/acled_attacks_district_year.csv`, is included) and `python/11_review3_acled_attacks.py`. The
authors hold the file; its checksum is in `environment/`.

## Software

All results come from Python 3.13 (statsmodels 0.15; versions in `environment/python_packages.txt`). The `R/` folder
in the authors' working directory is the superseded version 0.1 pipeline and is not part of this package.

## Pre-specification

`review3/PRESPECIFICATION_review3_checks.md` (third review), `review4/CORRECTION_zero_dose_definition.md` (20:17:01
UTC, written before the corrected rerun) and `review4/ADDENDUM_sensitivity_definitions.md` (20:59:54 UTC), each with
timestamp and SHA-256.
