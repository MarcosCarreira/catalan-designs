"""fz5_padic_check.py -- checks of the local estimates of A. Fauzan, "zeta(5) is irrational" (Zenodo, 17 Sep 2026), on the exact
determinants at n = 1, 2, 3 (K = 40n, N = 3n, h = 37n).

Run:  python3 fz5_padic_check.py        (no arguments; every switch is set below)
Needs: python-flint, mpmath, sympy; fz5_check.py in the same folder (its construction is imported).
Writes fz5_padic_check_log_<stamp>.txt and fz5_padic_check_out_<stamp>.json next to this file.
Expected runtime with 2 workers: about 5 minutes.

For every prime p < 3K, with Delta_K the exact determinant and F_K = S_K Delta_K:
  G1  (3.12): v_p^G(F_K) >= -6h floor(log_p(5K)) - h v_p(24);
  G2  Proposition 4.3, second part: for p > K, v_p^G(Delta_K) >= 0;
  G3  Proposition 4.3, first part: for p satisfying (4.9) (p >= 7, p <= K < 3p, p^2 > 2K, 2N < p, 5N <= 2p - 2),
      v_p^G(Delta_K) >= gamma_p^out of (4.14).
v_p^G is the minimum over the coefficients of the polynomial in X. The log lists, for every prime, the measured valuation and
the bound, so the sharpness of each bound can be read off.
"""
import os, sys, math, time, json, datetime
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE); sys.path.insert(0, HERE)
from multiprocessing import Pool
from flint import fmpq, fmpq_poly, fmpz_mat
from sympy import primerange

# ---------------- switches ----------------
WORKERS = 2
NS = [1, 2, 3]
# ------------------------------------------

STAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
LOGF = os.path.join(HERE, f"fz5_padic_check_log_{STAMP}.txt"); OUTF = os.path.join(HERE, f"fz5_padic_check_out_{STAMP}.json")
T0 = time.time(); _log = None
def say(s):
    global _log
    line = f"[{time.time()-T0:7.1f}s] {s}"; print(line, flush=True)
    if _log is None: _log = open(LOGF, "w")
    _log.write(line + "\n"); _log.flush()
FAILS = {}
def gate(name, ok, info=""):
    FAILS[name] = FAILS.get(name, 0) + (0 if ok else 1)
    if not ok: say(f"FAIL {name} {info}")

_src = open(os.path.join(HERE, 'fz5_check.py')).read().split('if __name__')[0]
_g = {'__file__': os.path.join(HERE, 'fz5_check.py'), '__name__': 'fz5mod'}
_saved = sys.argv; sys.argv = ['x']; exec(_src, _g); sys.argv = _saved
entries = _g['entries']

def vp_int(x, p):
    x = abs(int(x)); v = 0
    if x == 0: return 10**9
    while x % p == 0: x //= p; v += 1
    return v
def vq(q, p): return vp_int(q.p, p) - vp_int(q.q, p)
def vG(P, p): return min(vq(P[k], p) for k in range(P.degree()+1) if P[k] != 0)

_AI = _BI = None; _H = 0
def _init(ai, bi, h):
    global _AI, _BI, _H
    _AI, _BI, _H = ai, bi, h
def _det(x):
    h = _H
    return x, int(fmpz_mat(h, h, [_AI[i+j] + x*_BI[i+j] for i in range(h) for j in range(h)]).det())

def run(n):
    K, N, h, a, b, lead = entries(n)
    L = 1
    for q in a + b: L = L*int(q.q)//math.gcd(L, int(q.q))
    ai = [int((q*L).p) for q in a]; bi = [int((q*L).p) for q in b]
    xs = list(range(h+1)); ys = {}
    with Pool(WORKERS, initializer=_init, initargs=(ai, bi, h)) as pool:
        for x, d in pool.imap_unordered(_det, xs): ys[x] = d
    coef = [fmpq(ys[x]) for x in xs]
    for k in range(1, h+1):
        for i in range(h, k-1, -1): coef[i] = (coef[i]-coef[i-1]) / (xs[i]-xs[i-k])
    t = fmpq_poly([0, 1]); P = fmpq_poly([coef[h]])
    for k in range(h-1, -1, -1): P = P*(t - xs[k]) + coef[k]
    Delta = P / fmpq(L)**h
    fac = math.factorial
    SK = fmpq(fac(K)**(2*h) * 4**(h-1), fac(N)**(12*h) * math.prod(fac(2*i)**2 for i in range(1, h)))
    F = Delta * SK
    say(f"n = {n}: K = {K}, N = {N}, h = {h}; exact Delta_K of degree {Delta.degree()}")
    rows = []
    for p in primerange(2, 3*K):
        vF = vG(F, p); vD = vG(Delta, p)
        b312 = -6*h*int(math.floor(math.log(5*K)/math.log(p) + 1e-12)) - h*vp_int(24, p)
        gate("G1", vF >= b312, f"n={n} p={p} vG(F)={vF} bound {b312}")
        row = dict(p=p, vG_F=vF, bound_312=b312, vG_Delta=vD)
        txt = f"   p={p:4d}  vG(F_K)={vF:7d} >= {b312:7d}"
        if p > K:
            gate("G2", vD >= 0, f"n={n} p={p} vG(Delta)={vD}"); txt += f"  | p > K: vG(Delta_K)={vD} >= 0"
        elif p >= 7 and p <= K < 3*p and p*p > 2*K and 2*N < p and 5*N <= 2*p - 2:
            v = K - p*(K//p); u = max(0, N+v-p+1); tp = min(N, v) + u; rp = max(0, K+4*N-2*p+2)
            go = (-7*(K-p)+6*tp-1-min(rp, p-1-N+u)) if K < 2*p else (-7*(K-p)+3+12*N+5*tp-min(rp, p+u))
            gate("G3", vD >= go, f"n={n} p={p} vG(Delta)={vD} gamma_out={go}")
            row['gamma_out'] = go; txt += f"  | (4.9): vG(Delta_K)={vD} >= gamma_out={go} (excess {vD-go})"
        rows.append(row); say(txt)
    return dict(n=n, K=K, N=N, h=h, primes=rows)

if __name__ == "__main__":
    say(f"fz5_padic_check.py  workers={WORKERS}  n in {NS}")
    out = [run(n) for n in NS]
    for r in out:
        ex = [x['vG_Delta'] - x['gamma_out'] for x in r['primes'] if 'gamma_out' in x]
        say(f"n = {r['n']}: {len(ex)} primes satisfy (4.9); excess of v_p(Delta_K) over gamma_out: min {min(ex)}, max {max(ex)}")
    total = sum(FAILS.values())
    json.dump(dict(stamp=STAMP, gates=FAILS, results=out, gate_fails=total), open(OUTF, "w"), indent=1)
    say(f"gates: {FAILS}")
    say(f"GATE_FAILS = {total}")
    say(f"log saved to {LOGF}")
    say(f"json saved to {OUTF}")
