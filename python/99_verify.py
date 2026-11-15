"""99  Check every number printed in the manuscript against the outputs of the
current run.

manuscript_numbers.csv lists each number as it is printed, and where. A number
passes if the recomputed value rounds to it (within half a unit of its last
printed decimal). After any change to data or code, rerun everything and read
this report: every FAIL is a place in the manuscript that needs editing.
Exit code 1 if anything fails.
"""
import glob, os, sys
import pandas as pd
from common import KEYS

here = os.path.dirname(os.path.abspath(__file__))
exp = pd.read_csv(os.path.join(here, "manuscript_numbers.csv"), dtype=str)
got = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{KEYS}/*.csv"))])
got = dict(zip(got.key, got.value))
fast = got.get("ri_B", 20000) < 20000

ok, bad, miss, skip = [], [], [], []
for r in exp.itertuples():
    if r.key not in got:
        miss.append(r)
        continue
    if r.key == "ri_p" and fast:
        skip.append(r)
        continue
    v, s = got[r.key], r.manuscript_value
    dec = len(s.split(".")[1]) if "." in s else 0
    tol = 0.5 * 10 ** -dec + 1e-9
    passed = (v <= float(s) + 1e-9) if r.rule == "le" else abs(v - float(s)) <= tol
    (ok if passed else bad).append((r, v))

print(f"Checked {len(exp)} numbers printed in the manuscript against this run.\n")
for r, v in bad:
    print(f"  FAIL  {r.key:<40} manuscript {r.manuscript_value:<8} now {v:.4f}   "
          f"-> {r.where_in_manuscript}")
for r in miss:
    print(f"  MISSING  {r.key:<37} (no script produced it)  -> {r.where_in_manuscript}")
for r in skip:
    print(f"  SKIPPED  {r.key:<37} FAST=1 run uses 1,000 permutations; rerun without FAST")
print(f"\n{len(ok)} pass, {len(bad)} fail, {len(miss)} missing, {len(skip)} skipped")
sys.exit(1 if bad or miss else 0)
