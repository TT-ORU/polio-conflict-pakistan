"""Compare a fresh run (outputs/) with the frozen results (frozen_results/): every CSV in tables/ and keys/ must match
to 1e-9 and every figure PNG pixel for pixel. Run after run_all.sh."""
import glob, os, sys
import numpy as np, pandas as pd
from PIL import Image
os.chdir(os.path.dirname(os.path.abspath(__file__)))
ok, bad = 0, []
def same_csv(a, b):
    x, y = pd.read_csv(a), pd.read_csv(b)
    if list(x.columns) != list(y.columns) or x.shape != y.shape: return False
    for c in x.columns:
        if pd.api.types.is_numeric_dtype(x[c]) and pd.api.types.is_numeric_dtype(y[c]):
            if not np.allclose(x[c], y[c], rtol=1e-9, atol=1e-12, equal_nan=True): return False
        elif not (x[c].fillna("").astype(str) == y[c].fillna("").astype(str)).all(): return False
    return True
for f in sorted(glob.glob("frozen_results/tables/*.csv") + glob.glob("frozen_results/keys/*.csv")):
    g = f.replace("frozen_results/", "outputs/"); r = os.path.exists(g) and same_csv(f, g); ok += r
    if not r: bad.append(g)
for f in sorted(glob.glob("frozen_results/figures/*.png")):
    g = f.replace("frozen_results/", "outputs/")
    r = os.path.exists(g) and np.array_equal(np.asarray(Image.open(f)), np.asarray(Image.open(g))); ok += r
    if not r: bad.append(g)
print(f"{ok} files identical to the frozen version; {len(bad)} differ or are missing")
for b in bad: print("  differs:", b)
sys.exit(1 if bad else 0)
