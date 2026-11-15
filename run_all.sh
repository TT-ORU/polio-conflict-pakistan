#!/usr/bin/env bash
# Reproduce every analysis in the manuscript, then check every printed number.
#   bash run_all.sh            full run, about 10 minutes
#   FAST=1 bash run_all.sh     about 1 minute; randomisation inference uses
#                              1,000 permutations and its p-value is skipped
set -euo pipefail
# Two steps need the raw ACLED event file, which is not distributed (ACLED terms):
#   python3 python/00a_build_attack_counts.py /path/to/acled_pakistan_events_notes.csv
#   python3 python/11_review3_acled_attacks.py /path/to/acled_pakistan_events_notes.csv
# 00a writes data/acled_attacks_district_year.csv, which is included.
cd "$(dirname "$0")/python"
if [ ! -f ../data/pakistan_districts_gadm_level3.gpkg ]; then
  echo "District boundaries missing. Run first:  python3 python/get_boundaries.py"
  echo "(GADM boundaries may not be redistributed, so they are not in this repository.)"
  exit 1
fi
mkdir -p ../outputs/logs
for s in 00_build_panel 01_main_models 02_robustness 03_randomisation_inference \
         04_hurdle_and_temporal_order 05_temporal_validation 06_predicted_rates \
         07_spatial 08_figures 10_review3_checks 12_table1 13_figure1 14_figure2; do
  echo ">>> $s"
  python3 -W ignore "$s.py" 2>&1 | tee "../outputs/logs/$s.log"
done
echo ">>> 99_verify"
python3 -W ignore 99_verify.py 2>&1 | tee ../outputs/logs/99_verify.log
