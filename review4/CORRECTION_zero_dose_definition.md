# Correction of the district zero-dose rate (manuscript version 0.4)

Written on 21 October 2026, before the pipeline was rerun with the corrected input.

**Problem.** The district zero-dose rate used in versions up to 0.3 (data/district_zero_dose_mics_v03_original.csv)
did not count children whose mother reported that they had never received any vaccination (MICS6 IM11 = 2). For
these children the vaccine-specific recall questions are skipped by design, so they were treated as undetermined
rather than zero-dose. This under-counted zero-dose children, most in the districts with the lowest coverage. The
error was found while building the companion zero-dose paper.

**Correction.** data/district_zero_dose_mics.csv is replaced by the district rate from the companion paper
(12_zero_dose/18_polio_zero_dose.py): weighted proportion of children aged 12-23 months with no pentavalent dose by
card or recall, counting never-vaccinated reports as zero-dose (when such a report conflicts with a positive
pentavalent recall answer, the specific answer is used). MICS districts are mapped to panel districts through the
same crosswalk; a panel district linked to several MICS districts pools their children with survey weights. The
age group is 12-23 months, the standard zero-dose definition, so both papers use the same measure.

**What is rerun.** Every script in run_all.sh, unchanged apart from the input file and the frozen-panel check, which
now expects zero_dose_rate_mics to differ. All pre-specified analyses, the review-3 checks and the temporal
validation are rerun as before; no new analyses are added. Tertile cut-points, flagged districts and every number
that depends on the zero-dose rate are expected to change, and 99_verify.py lists each manuscript number to update.
Results are reported whatever they show.

Old and new rates: 120 districts; correlation 0.866; mean 0.210 (old) and 0.268 (new).
