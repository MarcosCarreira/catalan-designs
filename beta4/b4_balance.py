"""b4_balance.py -- measured balance of the H2 half-integer determinants, weight 4 (beta(4)) against
weight 2 (Catalan's constant G, calibration), for a few designs and n = 6 .. 40.

Run:  python3 b4_balance.py          (no arguments; every switch is set below)
Writes b4_balance_log_<stamp>.txt and b4_balance_out_<stamp>.json next to this file.

Construction (as in h2_scout.py, Section 6 of the extension paper):
  Z'_m(i,j) = sum_{u>=0} (-1)^u / ((i+u+1/2)^m (j+u+1/2)),  i, j even.
  Rows a3_r = coefficients of t^{2(alpha+r)} (1-t^2)^h, columns c_k = coefficients of s^{2(b+k)} (1-s^2)^q.
  W=4: a2, a1 on the column support chosen so that G and pi^3 cancel; entries in Q + Q beta(4).
  W=2: the single array a3 with m=1; entries in Q + Q G.
  Design ratios (ra, rh, rb, rq) give alpha, h, b, q = round(ratio * n).
Measured per (design, W, n):
  L_true = log|Delta| / n^2                         (Delta computed with the true constant)
  L_low  = -log den(Delta_hyp) / n^2                (exact Delta with the constant replaced by 1/3)
  gap    = L_true - L_low                            (an irrationality proof in this family needs gap < 0)
Gates (pre-registered; GATE_FAILS counts failures):
  P1  positive control: decomposition = directly summed series, |diff| < 1e-30, W=4 and W=2
  P2  unwanted-constant residues = 0 in every run
  P3  Delta (true) != 0 and Delta_hyp != 0 in every run
  P4  the mpmath determinant of the hypothetical matrix agrees with the exact Fraction value,
      relative error < 1e-30 (checks that dps is enough for the cancellation in the determinant)
Read-out (printed at the end, no gate): per design and weight, the gap sequence, a fit
  gap(n) = c0 + c1/n + c2/n^2 over n >= 12, and the ratio c0(W=4)/c0(W=2).
Expected runtime on the Mac with 2 workers: about 4-8 minutes (n=40 dominates, about 1.5 min per W=4 case).
"""
import os, sys, math, time, json, datetime
os.environ.setdefault("OMP_NUM_THREADS", "2"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
from fractions import Fraction as F
from math import comb
from multiprocessing import Pool
import mpmath as mp

# ---------------- switches ----------------
WORKERS = 2
DESIGNS = [  # (ra, rh, rb, rq)
    (0.0, 0.0,   0.0, 0.0),
    (0.0, 0.125, 0.0, 0.125),
    (0.0, 0.25,  0.0, 0.25),
    (0.0, 0.0,   0.0, 0.25),
    (0.25, 0.25, 0.0, 0.25),
]
WEIGHTS = [4, 2]
NS = [6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40]
FIT_NMIN = 12
P1_TOL = mp.mpf("1e-30")
P4_TOL = mp.mpf("1e-30")
# ------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))
STAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
LOGF = os.path.join(HERE, f"b4_balance_log_{STAMP}.txt")
OUTF = os.path.join(HERE, f"b4_balance_out_{STAMP}.json")
T0 = time.time()
_log = None
def say(s):
    global _log
    line = f"[{time.time()-T0:7.1f}s] {s}"
    print(line, flush=True)
    if _log is None: _log = open(LOGF, "w")
    _log.write(line + "\n"); _log.flush()

# ---------- exact decomposition ----------
_R = {}
def R(j, m):  # sum_{v<j} (-1)^v / (v+1/2)^m
    if (j, m) not in _R:
        _R[(j, m)] = sum((F((-1)**v) / F(2*v+1, 2)**m for v in range(j)), F(0))
    return _R[(j, m)]
_Z = {}
def Zrat(m, i, j):
    """Z'_m(i,j) = out[0] + sum_k out[k] beta(k), i and j even (beta(1) = pi/4, beta(2) = G, beta(3) = pi^3/32)."""
    key = (m, i, j)
    if key in _Z: return _Z[key]
    out = [F(0)]*5
    def addT(mm, idx, c):  # c * T_mm(idx), T_mm(idx) = 2^mm beta(mm) - R(idx, mm) for even idx
        out[mm] += c * 2**mm; out[0] -= c * R(idx, mm)
    if i == j:
        addT(m+1, i, F(1))
    else:
        d = i - j
        addT(1, j, F(1, d**m)); addT(1, i, F(-1, d**m))
        if m >= 2: addT(2, i, F(-1, d**(m-1)))
        if m >= 3: addT(3, i, F(-1, d**(m-2)))
    _Z[key] = out
    return out

def arrays_for_row(a3, J, W):
    if W == 2: return ((1, a3),)
    a2 = {j: sum((v / (i-j) for i, v in a3.items() if i != j), F(0)) for j in J}
    a1 = {j: sum((v / (i-j) for i, v in a2.items() if i != j), F(0))
             + sum((v / F(i-j)**2 for i, v in a3.items() if i != j), F(0)) for j in J}
    return ((3, a3), (2, a2), (1, a1))

def build(n, alpha, h, b, q, W):
    rows = [{2*(alpha+r+l): F((-1)**l * comb(h, l)) for l in range(h+1)} for r in range(n)]
    cols = [{2*(b+k+l): F((-1)**l * comb(q, l)) for l in range(q+1)} for k in range(n)]
    J = sorted(set(j for c in cols for j in c))
    tgt = 4 if W == 4 else 2
    Frat = [[F(0)]*n for _ in range(n)]; Ftg = [[F(0)]*n for _ in range(n)]; resid = 0
    for r, a3 in enumerate(rows):
        arrs = arrays_for_row(a3, J, W)
        for k, c in enumerate(cols):
            tot = [F(0)]*5
            for j, cv in c.items():
                for m, arr in arrs:
                    for i, av in arr.items():
                        if av == 0: continue
                        z = Zrat(m, i, j); f = av*cv
                        for key in range(5):
                            if z[key]: tot[key] += f*z[key]
            if any(tot[x] != 0 for x in (1, 2, 3, 4) if x != tgt): resid += 1
            Frat[r][k] = tot[0]; Ftg[r][k] = tot[tgt]
    return Frat, Ftg, resid

def det_exact(M):
    """exact determinant of a Fraction matrix: scale rows to integers, Bareiss."""
    n = len(M); scale = 1; A = []
    for row in M:
        L = 1
        for x in row: L = L * x.denominator // math.gcd(L, x.denominator)
        scale *= L; A.append([int(x * L) for x in row])
    sign = 1; prev = 1
    for c in range(n-1):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None: return F(0)
        if p != c: A[c], A[p] = A[p], A[c]; sign = -sign
        for r in range(c+1, n):
            for k in range(c+1, n):
                A[r][k] = (A[r][k]*A[c][c] - A[r][c]*A[c][k]) // prev
            A[r][c] = 0
        prev = A[c][c]
    return F(sign * A[n-1][n-1], scale)

def mpq(x): return mp.mpf(x.numerator) / x.denominator

def const_of(W):
    return (mp.zeta(4, mp.mpf(1)/4) - mp.zeta(4, mp.mpf(3)/4)) / 256 if W == 4 else +mp.catalan

def run(task):
    des, W, n = task
    t = time.time()
    ra, rh, rb, rq = des
    alpha, h, b, q = round(ra*n), round(rh*n), round(rb*n), round(rq*n)
    Frat, Ftg, resid = build(n, alpha, h, b, q, W)
    mp.mp.dps = 60 + 4*n*n
    c = const_of(W)
    M = mp.matrix(n, n); Mh = mp.matrix(n, n)
    third = mp.mpf(1)/3
    for r in range(n):
        for k in range(n):
            x, y = mpq(Frat[r][k]), mpq(Ftg[r][k])
            M[r, k] = x + c*y; Mh[r, k] = x + third*y
    dtrue = mp.det(M)
    Dh = det_exact([[Frat[r][k] + F(1, 3)*Ftg[r][k] for k in range(n)] for r in range(n)])
    dh_mp = mp.det(Mh)
    p4 = (abs(dh_mp - mpq(Dh)) / abs(mpq(Dh))) if Dh != 0 else mp.inf
    Lt = float(mp.log(abs(dtrue)) / n**2) if dtrue != 0 else None
    Ll = float(-mp.log(Dh.denominator) / n**2) if Dh != 0 else None
    return dict(design=list(des), W=W, n=n, alpha=alpha, h=h, b=b, q=q, resid=resid,
                det_true_nonzero=bool(dtrue != 0), det_hyp_nonzero=bool(Dh != 0),
                p4_relerr=float(p4) if p4 != mp.inf else None, p4_ok=bool(p4 < P4_TOL),
                L_true=Lt, L_low=Ll, gap=(Lt - Ll) if (Lt is not None and Ll is not None) else None,
                secs=time.time() - t)

def positive_control(W):
    """decomposition vs directly summed series on four entries of a small design."""
    mp.mp.dps = 50
    c = const_of(W)
    n, alpha, h, b, q = 3, 1, 2, 0, 2
    Frat, Ftg, resid = build(n, alpha, h, b, q, W)
    rows = [{2*(alpha+r+l): F((-1)**l*comb(h, l)) for l in range(h+1)} for r in range(n)]
    cols = [{2*(b+k+l): F((-1)**l*comb(q, l)) for l in range(q+1)} for k in range(n)]
    J = sorted(set(j for cc in cols for j in cc))
    worst = mp.mpf(0)
    for r in (0, 2):
        arrs = arrays_for_row(rows[r], J, W)
        for k in (0, 2):
            def term(u):
                s = mp.mpf(0)
                for j, cv in cols[k].items():
                    for m, arr in arrs:
                        for i, av in arr.items():
                            s += mpq(av)*cv / ((i+u+mp.mpf(1)/2)**m * (j+u+mp.mpf(1)/2))
                return (-1)**int(u) * s
            direct = mp.nsum(term, [0, mp.inf])
            dec = mpq(Frat[r][k]) + c*mpq(Ftg[r][k])
            worst = max(worst, abs(direct - dec))
    return worst, resid

def fit(ns, gs):
    """least squares gap = c0 + c1/n + c2/n^2"""
    import numpy as np
    A = np.array([[1, 1/n, 1/n**2] for n in ns]); y = np.array(gs)
    sol = np.linalg.lstsq(A, y, rcond=None)[0]
    return [float(v) for v in sol]

if __name__ == "__main__":
    fails = 0
    say(f"b4_balance.py  designs={DESIGNS}  W={WEIGHTS}  n={NS}  workers={WORKERS}")
    for W in WEIGHTS:
        worst, resid = positive_control(W)
        ok = worst < P1_TOL and resid == 0
        fails += (not ok)
        say(f"P1 W={W}: max |direct - decomposition| = {mp.nstr(worst, 3)}, residues {resid}  -> {'ok' if ok else 'FAIL'}")
    tasks = [(d, W, n) for n in sorted(NS, reverse=True) for d in DESIGNS for W in WEIGHTS]
    results = []
    with Pool(WORKERS) as pool:
        for res in pool.imap_unordered(run, tasks):
            results.append(res)
            g = res["gap"]
            flag = []
            if res["resid"] != 0: flag.append("P2 FAIL")
            if not (res["det_true_nonzero"] and res["det_hyp_nonzero"]): flag.append("P3 FAIL")
            if not res["p4_ok"]: flag.append("P4 FAIL")
            say(f"design={tuple(res['design'])} W={res['W']} n={res['n']:2d} (alpha,h,b,q)=({res['alpha']},{res['h']},{res['b']},{res['q']}): "
                f"L_true={res['L_true']:.5f} L_low={res['L_low']:.5f} gap={g:.5f}  P4 relerr={res['p4_relerr']:.1e}  "
                f"({res['secs']:.1f}s) {' '.join(flag)}" if g is not None else
                f"design={tuple(res['design'])} W={res['W']} n={res['n']}: determinant zero {' '.join(flag)}")
    fails += sum(1 for r in results if r["resid"] != 0)
    fails += sum(1 for r in results if not (r["det_true_nonzero"] and r["det_hyp_nonzero"]))
    fails += sum(1 for r in results if not r["p4_ok"])
    say("---- read-out (no gate) ----")
    summary = []
    for d in DESIGNS:
        row = {"design": list(d)}
        for W in WEIGHTS:
            rs = sorted((r for r in results if tuple(r["design"]) == d and r["W"] == W and r["gap"] is not None),
                        key=lambda r: r["n"])
            seq = [(r["n"], r["gap"]) for r in rs]
            pts = [(n, g) for n, g in seq if n >= FIT_NMIN]
            c = fit([p[0] for p in pts], [p[1] for p in pts]) if len(pts) >= 3 else None
            row[f"W{W}"] = {"gaps": seq, "fit": c}
            say(f"design={d} W={W}: gaps " + ", ".join(f"n={n}:{g:.3f}" for n, g in seq)
                + (f"  | fit c0={c[0]:.3f} c1={c[1]:.2f} c2={c[2]:.1f}" if c else ""))
        f4, f2 = row["W4"]["fit"], row["W2"]["fit"]
        if f4 and f2 and f2[0] != 0:
            row["ratio_c0"] = f4[0]/f2[0]
            say(f"design={d}: limit gap W=4 {f4[0]:.3f}, W=2 {f2[0]:.3f}, ratio {f4[0]/f2[0]:.2f}")
        summary.append(row)
    json.dump({"stamp": STAMP, "designs": DESIGNS, "weights": WEIGHTS, "ns": NS,
               "results": results, "summary": summary, "gate_fails": fails}, open(OUTF, "w"), indent=1)
    say(f"GATE_FAILS = {fails}")
    say(f"log saved to {LOGF}")
    say(f"json saved to {OUTF}")
