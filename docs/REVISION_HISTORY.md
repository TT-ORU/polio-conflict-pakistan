
## Fourth revision (version 0.4): corrected zero-dose rate

- **What changed.** The district zero-dose rate used up to version 0.3 omitted children whose mother reported that
  they had never received any vaccination (MICS6 IM11 = 2; the vaccine-specific questions are skipped for them). The
  corrected rate counts them as zero-dose and uses children aged 12–23 months, the definition of the companion
  zero-dose paper (12_zero_dose/22_polio_zero_dose_corrected.py). `data/district_zero_dose_mics.csv` now holds the
  corrected rate; the previous one is kept as `data/district_zero_dose_mics_v03_original.csv`, and
  `data/zero_dose_check_review4.csv` holds corrected 12–23, 12–35 and 0–35 month rates and the previous rate.
- **Specification.** `review4/CORRECTION_zero_dose_definition.md` (20:17:01 UTC, sha256 in `review4/CORRECTION.sha256`)
  was written before the rerun; `review4/ADDENDUM_sensitivity_definitions.md` (20:59:54 UTC) sets the Table S10 rows.
- **Code.** `00_build_panel.py` now expects `zero_dose_rate_mics` to differ from the frozen panel; `10_review3_checks.py`
  compares the corrected and previous definitions; `08_figures.py` computes the Figure 4 note from the outputs.
- **Results.** Full rerun (`run_all_console_v4.log`, 20,000 permutations); `99_verify.py`: 367 of 367 manuscript
  numbers pass against `manuscript/Polio_conflict_manuscript_v0.4_zero_dose_TRACKED.docx`. Main changes: zero-dose IRR
  1.83 to 1.70; interaction 1.098 to 1.133 (randomisation-inference p 0.046 to 0.015); the high/high tertile cell has 21
  districts (was 20), without Karachi West; temporal-validation IRR 7.11 to 9.41; stratified estimate in the top
  zero-dose tertile no longer excludes 1; an inverse violence association appears where zero-dose prevalence is low.

### Software and the R folder (version 0.6)

Every number in the manuscript comes from the Python pipeline in `python/` (run with `run_all.sh`; versions in
`requirements.txt`) and is checked by `python/99_verify.py`; `review4/verify_tables_v05.py` checks every table cell of
the manuscript against `outputs/`. The `R/` folder holds the original version 0.1 pipeline (cases to 2024, event-study
specification) and is not used for any reported result; it is kept only as a record.

### Prior-history adjustment in the temporal validation

Table S8 reports two adjustments of the high-risk flag for a district's own 2018–2021 history: log(1 + prior cases),
IRR 3.36 (1.44–7.84), and log(1 + prior cases per million person-years), IRR 2.72 (1.09–6.81). The abstract reports the
second (prior incidence) and gives the first alongside it.
