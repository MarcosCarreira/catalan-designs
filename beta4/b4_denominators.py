"""b4_denominators.py -- where the extra denominator of the weight-4 half-integer determinants comes from.

Run:  python3 b4_denominators.py     (no arguments; needs b4_balance.py in the same folder)
Writes b4_denominators_log_<stamp>.txt and b4_denominators_out_<stamp>.json next to this file.

For each design and n, three exact determinants, every constant (pi/4, G, pi^3/32, beta(4)) replaced by 1/3:
  W2      weight 2: array Z'_1 with the rows a3                         (the Catalan calibration of b4_balance.py)
  W4-a3   weight 4: array Z'_3 with the rows a3 only (G and pi^3 NOT cancelled)
  W4      weight 4: arrays Z'_3, Z'_2, Z'_1 with a3, a2, a1 (G and pi^3 cancelled, as in b4_balance.py)
Printed: -log den(Delta)/n^2 for each. If W4-a3 and W4 agree, the cancellation of G and pi^3 adds no denominator,
and the difference to W2 comes from the weight.
Gates: D1 residues of unwanted constants = 0 for W2 and W4 (W4-a3 keeps them by construction); D2 all determinants
nonzero.  Expected runtime with 2 workers: about 1 minute.
"""
import os, sys, math, time, json, datetime
os.environ.setdefault("OMP_NUM_THREADS", "2"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
from fractions import Fraction as F
from math import comb
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from b4_balance import Zrat, arrays_for_row, det_exact

# ---------------- switches ----------------
WORKERS = 2
DESIGNS = [(0.0, 0.0, 0.0, 0.0), (0.0, 0.25, 0.0, 0.25)]   # (ra, rh, rb, rq); ra = rb = 0 as in b4_balance.py
NS = [10, 14, 20, 28, 40]
VARIANTS = ["W2", "W4-a3", "W4"]
# ------------------------------------------

STAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
LOGF = os.path.join(HERE, f"b4_denominators_log_{STAMP}.txt")
OUTF = os.path.join(HERE, f"b4_denominators_out_{STAMP}.json")
T0 = time.time(); _log = None
def say(s):
    global _log
    line = f"[{time.time()-T0:7.1f}s] {s}"; print(line, flush=True)
    if _log is None: _log = open(LOGF, "w")
    _log.write(line + "\n"); _log.flush()

def run(task):
    des, var, n = task
    t = time.time()
    ra, rh, rb, rq = des
    alpha, h, b, q = round(ra*n), round(rh*n), round(rb*n), round(rq*n)
    rows = [{2*(alpha+r+l): F((-1)**l*comb(h, l)) for l in range(h+1)} for r in range(n)]
    cols = [{2*(b+k+l): F((-1)**l*comb(q, l)) for l in range(q+1)} for k in range(n)]
    J = sorted(set(j for c in cols for j in c))
    W = 2 if var == "W2" else 4
    tgt = 2 if W == 2 else 4
    M = []; resid = 0
    for a3 in rows:
        arrs = arrays_for_row(a3, J, W)
        if var == "W4-a3": arrs = arrs[:1]
        row = []
        for c in cols:
            tot = [F(0)]*5
            for j, cv in c.items():
                for m, arr in arrs:
                    for i, av in arr.items():
                        if av == 0: continue
                        z = Zrat(m, i, j); f = av*cv
                        for key in range(5):
                            if z[key]: tot[key] += f*z[key]
            if any(tot[x] != 0 for x in (1, 2, 3, 4) if x != tgt): resid += 1
            row.append(tot[0] + F(1, 3)*(tot[1] + tot[2] + tot[3] + tot[4]))
        M.append(row)
    D = det_exact(M)
    return dict(design=list(des), variant=var, n=n, resid=resid, nonzero=bool(D != 0),
                L_low=(-math.log(D.denominator)/n**2) if D != 0 else None, secs=time.time()-t)

if __name__ == "__main__":
    say(f"b4_denominators.py  designs={DESIGNS}  n={NS}  variants={VARIANTS}  workers={WORKERS}")
    tasks = [(d, v, n) for n in sorted(NS, reverse=True) for d in DESIGNS for v in VARIANTS]
    res = []
    with Pool(WORKERS) as pool:
        for r in pool.imap_unordered(run, tasks):
            res.append(r)
            say(f"design={tuple(r['design'])} {r['variant']:6s} n={r['n']:2d}: -log den/n^2 = "
                + (f"{r['L_low']:.4f}" if r['L_low'] is not None else "det 0")
                + f"  residues {r['resid']}  ({r['secs']:.1f}s)")
    fails = sum(1 for r in res if r["variant"] != "W4-a3" and r["resid"] != 0)
    fails += sum(1 for r in res if not r["nonzero"])
    say("---- table: -log den(Delta)/n^2 ----")
    for d in DESIGNS:
        for n in NS:
            v = {r["variant"]: r["L_low"] for r in res if tuple(r["design"]) == d and r["n"] == n}
            say(f"design={d} n={n:2d}: " + "  ".join(f"{k}={v[k]:.3f}" for k in VARIANTS if v.get(k) is not None))
    json.dump({"stamp": STAMP, "designs": DESIGNS, "ns": NS, "results": res, "gate_fails": fails}, open(OUTF, "w"), indent=1)
    say(f"GATE_FAILS = {fails}")
    say(f"log saved to {LOGF}")
    say(f"json saved to {OUTF}")
