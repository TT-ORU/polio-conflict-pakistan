# Pre-specification: checks requested in the third review of the polio paper

Written on 19 October 2026, before any of the models below was fitted. At the time of writing the
published results (Tables 2 to 4, S1 to S9) were known; none of the checks below had been run.
Saved before the analysis script is executed.

## A. Input checks (no models)

1. **Zero-dose variable.** Confirm that one definition (never received any Pentavalent dose, card or
   recall, MICS6 child questionnaire) is applied in all four provinces. Report the number of districts with
   a zero-dose value and which districts lack one.
2. **ACLED health-worker attacks.** Rebuild the event flags from the event notes, compare with the
   counts used in the panel, and read a random sample of 40 flagged events (seed 20261019) to estimate
   the share that are attacks on health or polio workers.
3. **Case reconstruction.** Compare the annual national totals with published GPEI totals.
4. **Population.** Check that GHSL district population changes from year to year and how years without
   a GHSL epoch are filled.
5. **Moran's I, 2021.** Inspect the permutation distribution for 2021, the year with one case.

## B. Analyses

All use the primary interaction model of Table 3 (Poisson, log population offset, province and year
fixed effects, district-clustered SEs, 840 district-years, 120 districts) unless stated.

1. **Included versus excluded districts.** Compare the 120 districts in the interaction sample with the
   21 outside it: cases, incidence, violence, population, province. Descriptive only.
2. **Targeting rules (temporal validation).** Using 2018 to 2021 data only, compare rules for flagging
   districts, outcome = reported cases 2022 to 2025 in the same 120 districts:
   (a) top violence tertile; (b) top zero-dose tertile; (c) both (the paper's rule);
   (d) any reported case 2018 to 2021; (e) the 19 districts with the most cases 2018 to 2021;
   (f) the 19 districts with the highest incidence 2018 to 2021; (g) (c) or (d).
   Report districts flagged, population share, share of later cases, and cases per million in flagged
   versus other districts. No significance testing; the comparison is descriptive.
3. **Interaction inference.**
   (a) Wild cluster score bootstrap (Kline and Santos), null imposed, Rademacher and Webb weights,
   9,999 draws, clustered by district.
   (b) Model refitted excluding 2024, and excluding 2024 and 2025.
   (c) Randomisation inference is not rerun; its null hypothesis is described in the response.
4. **Additional covariates.** Three district-level covariates, entered together and then one at a
   time, each interacted with nothing: travel time to the nearest city (2015), log(1 + OSM health
   facilities per 1,000 km²), and log night-time lights (VIIRS, mean 2019 to 2024). Question: does the
   interaction estimate change? The interaction is judged robust if its point estimate stays above 1.05
   with all three covariates included.
5. **Spatial inference.** Conley standard errors for the interaction with a uniform kernel at 100, 200
   and 300 km between district centroids, allowing correlation across all years within the cutoff.

Nothing above changes the primary analysis. All results are reported whatever they show.
