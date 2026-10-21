"""
cat_mu_search.py -- design search with the irrationality exponent as the objective.

OBJECTIVE (quick evaluator, numerical, not certified): for a design d = (n0, b0, q0, g0, h0)
    gamma^  = f_low(d) - U_2^(d)                (U_2^ = kappa = 2 value of solve(), 300 coefficients),
    Sigma_2^ = sup X + sup Y of that solution (grid suprema),
    rho_1^  = (2 - alpha - gg) log 2 - Sigma_2^ + 2 log 2,
    rho_2^  = rho_1^ + (eta + gq) log 2 + (1 - gg) log(1 + sqrt 2),
    mu^     = rho_2^ ((1 - gg)/gamma^ + n0/2)  (Proposition xprop:power with Lemma xlem:nonzero for
                                                 p^2 not dividing b; paper v2 of 20 Oct 2026),
with alpha = (b0+q0)/n0, gg = g0/n0, gq = q0/n0, eta = h0/n0; mu^ = +inf if gamma^ <= 0.
f_low, solve() and the grid are those of cat_design_search.py (copied below, solve() also returns
Sigma_2^). The certified values come afterwards from run_design_cert.py for the chosen design.

SEARCH: for each n0 in N0_LIST, start at round(n0*(0.16, 0.06, 0.083, 0.053)) (ratios of the two known
good designs) and descend: evaluate the 8 neighbours (b0, q0, g0, h0 each -/+ 1), move to the best one
if it lowers mu^ by more than 1e-3, repeat until none does (at most MAXEV evaluations per n0).
Admissible: b0 >= 0, q0 >= 1, g0 >= 1, h0 >= 0, 2 g0 < n0 - 1 (C1), 2 b0 < H0 (C2).
The values of n0 are distributed over WORKERS processes (one BLAS thread each).

GATES (pre-registered)
  S1 the start design (48, 7, 3, 4, 2) of cat_design_search.py reproduces its score +0.001676
     within 5e-5                                                                       [0 fails]
  S2 control (96, 15, 6, 8, 5): gamma^ = +0.002443 within 5e-5 (cat_design_search.py record run) [0 fails]
REPORT: per n0 the descent path and its optimum; the overall ranking by mu^; mu^ of the controls
(96, 15, 6, 8, 5) and (36, 6, 2, 3, 2).

RUN: python3 cat_mu_search.py   (no arguments; needs numpy and scipy)
TIME: about 5 s per evaluation; 22 values of n0, about 15-30 evaluations each, 2 workers:
  expected 20-60 min.
OUTPUT: screen + cat_mu_search_log_<YYYYMMDD_HHMMSS>.txt + cat_mu_search_out.json, next to the script.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
import sys, json, time, math
from fractions import Fraction as Fr
import multiprocessing as mp
import numpy as np
from scipy.optimize import minimize

WORKERS = 2
N0_LIST = [12, 16, 18, 20, 24, 28, 30, 32, 36, 40, 42, 44, 48, 52, 54, 56, 60, 64, 72, 80, 88, 96]
MAXEV = 60
START_RATIOS = (0.16, 0.06, 0.083, 0.053)
CONTROLS = [(96, 15, 6, 8, 5), (36, 6, 2, 3, 2)]
LOG2 = math.log(2)
LOGS = math.log(1 + math.sqrt(2))

T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
STAMP = time.strftime('%Y%m%d_%H%M%S')
LOGPATH = os.path.join(HERE, 'cat_mu_search_log_' + STAMP + '.txt')
def log(s):
    line = f"[{time.time()-T0:7.1f}s] {s}"
    print(line, flush=True)
    with open(LOGPATH, 'a') as fh: fh.write(line + '\n')
FAILS = {}
def gate(name, ok, info=""):
    FAILS[name] = FAILS.get(name, 0) + (0 if ok else 1)
    if not ok: log(f"FAIL {name} {info}")

pos = lambda y: y if y > 0 else 0

def derived(n, b, q, g, h):
    a = b + q; L = n + a; C = L + g; H = C + h; A = a + 2*g
    return a, L, C, H, A

def Dfun(x, n, b, q, g, h):
    a, L, C, H, A = derived(n, b, q, g, h)
    if 2*x <= H:
        d0 = pos(min(x, L - x, A) - max(b, H - 2*x))
        return 2*n - pos(d0 - q)
    K, J = L - x, H - x
    R = pos(K) + pos(K - A) + pos(J - max(b, K))
    S = pos(min(A, K) - b) + pos(min(b, J) - max(0, K))
    return min(2*n, n + R, 2*R + S)

def breakpoints(n, b, q, g, h):
    """candidate kinks/jumps: all crossings among the linear forms that appear in D."""
    a, L, C, H, A = derived(n, b, q, g, h)
    consts = [0, b, q, A, n, 2*n, L, H, A + b, b + q, A + q, n - b, n - q]
    lins = [(1, 0), (1, q), (-1, L), (-1, L - q), (-2, H), (-2, H + q), (-1, H), (-1, L - A),
            (-1, H - b), (-1, L - b), (-2, H - b)]
    pts = set([Fr(0), Fr(H), Fr(H, 2)])
    forms = [(0, c) for c in consts] + lins
    for i in range(len(forms)):
        for j in range(i + 1, len(forms)):
            (s1, c1), (s2, c2) = forms[i], forms[j]
            if s1 != s2:
                x = Fr(c2 - c1, s1 - s2)
                if 0 <= x <= H: pts.add(x)
    return sorted(pts)

def integrate(n, b, q, g, h, exact=True):
    """integral of D over [0, H]. D is piecewise linear and may jump at isolated points.
    On each piece the interior is tested for linearity at 7 points; a linear piece contributes
    width * D(midpoint) (endpoint jumps do not matter). Pieces that are not linear are bisected;
    after 60 levels a piece is accepted with its midpoint value and counted as unresolved."""
    a, L, C, H, A = derived(n, b, q, g, h)
    num = (lambda v: Fr(v)) if exact else (lambda v: float(v))
    f = lambda x: num(Dfun(x, n, b, q, g, h))
    pts = [num(p) for p in breakpoints(n, b, q, g, h)]
    total = num(0); pieces = []; unresolved = 0
    stack = [(pts[i], pts[i+1], 0) for i in range(len(pts)-1)][::-1]
    while stack:
        lo, hi, depth = stack.pop()
        if hi <= lo: continue
        w = hi - lo
        xs = [lo + w*num(k)/8 for k in range(1, 8)]
        ys = [f(x) for x in xs]
        slope = (ys[6] - ys[0])/(xs[6] - xs[0])
        lin = all(abs(ys[k] - (ys[0] + slope*(xs[k]-xs[0]))) <= (0 if exact else 1e-9) for k in range(7))
        if not lin and depth < 60:
            mid = (lo + hi)/2; stack.append((mid, hi, depth+1)); stack.append((lo, mid, depth+1)); continue
        if not lin: unresolved += 1
        ymid = ys[3]
        total += w*ymid
        pieces.append((lo, hi, ys[0] - slope*(xs[0]-lo), ys[0] + slope*(hi-xs[0])))
    return total, pieces, unresolved

def f_low(n, b, q, g, h, exact=False):
    I, _, unres = integrate(n, b, q, g, h, exact=exact)
    gate("G1", unres == 0, f"{(n,b,q,g,h)}: {unres} unresolved pieces")
    delta = Fr(q + g + h, n)
    val = -float(I)/n**2 - (0.5 + float(delta)/2 + float(delta)**2/8)*LOG2
    return val, I


K = 300                     # number of trial coefficients per sequence
NX = 6000                   # optimisation grid
def design_ratios(n, b, q, g, h):
    a = b + q
    return a/n, b/n, q/n, g/n, h/n
th = np.linspace(0, np.pi, NX+2)[1:-1]
X = np.cos(th); X = X[np.abs(X) > 1e-6]; X = np.concatenate([[-1.0], X])
S = (1 - np.cos(th))/2
kk = np.arange(1, K+1)
def cheb(x):
    T = np.zeros((K+1, len(x))); T[0] = 1; T[1] = x
    for k in range(2, K+1): T[k] = 2*x*T[k-1] - T[k-2]
    return T[1:]
TX, TS = cheb(X), cheb(S)
PX = np.array([X**k for k in kk])
def setup(des, kap):
    al, be, gq, gg, et = design_ratios(*des)
    W = (al+2*gg)*np.log(np.abs(X)) + 2*et*np.log(1-X) - (kap/2+al+et+gg)*np.log(1+X**2)
    Dx = 2*gg - 2*X**2/(1+X**2)
    WS = be*np.log(S) + gq*np.log(1-S)
    const = (-1+al+gg)*np.log(2)
    return W, Dx, WS, const
def lse(z, tau):
    m = z.max(); w = np.exp((z-m)/tau); s = w.sum()
    return m + tau*np.log(s), w/s
def solve(des, kap, verbose=False):
    W, Dx, WS, const = setup(des, kap)
    nv = 2*K + (1 if kap == 1 else 0)
    def unpack(z):
        p, v = z[:K], z[K:2*K]; lam = z[2*K]**2 if kap == 1 else 0.0
        return p, v, lam
    def F(z, tau):
        p, v, lam = unpack(z)
        Xv = W + lam*Dx - 2*kap*(p/kk) @ TX - (v/kk) @ PX
        Yv = WS + 2*(v/kk) @ TS
        sx, wx = lse(Xv, tau); sy, wy = lse(Yv, tau)
        val = const + kap*np.sum(p*p/kk) + 0.5*np.sum(v*v/kk) + sx + sy
        gp = 2*kap*p/kk - 2*kap*(TX @ wx)/kk
        gv = v/kk - (PX @ wx)/kk + 2*(TS @ wy)/kk
        g = np.concatenate([gp, gv])
        if kap == 1: g = np.concatenate([g, [2*z[2*K]*(Dx @ wx)]])
        return val, g
    z = np.zeros(nv)
    if kap == 1: z[2*K] = 1.5   # lambda = 2.25 start (lambda = z^2 has zero gradient at z = 0)
    for tau in (1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5):
        r = minimize(F, z, args=(tau,), jac=True, method='L-BFGS-B', options={'maxiter': 3000})
        z = r.x
    p, v, lam = unpack(z)
    # exact-ish evaluation of the objective on a fine grid
    xf = np.cos(np.linspace(0, np.pi, 400002)[1:-1]); xf = xf[np.abs(xf) > 1e-7]; xf = np.concatenate([[-1.0], xf])
    sf = (1-np.cos(np.linspace(0, np.pi, 400002)[1:-1]))/2
    al, be, gq, gg, et = design_ratios(*des)
    def chebv(x, c):
        t0, t1 = np.ones_like(x), x; out = np.zeros_like(x)
        for k in range(1, K+1):
            out += c[k-1]*t1/k; t0, t1 = t1, 2*x*t1 - t0
        return out
    Xf = (al+2*gg)*np.log(np.abs(xf)) + 2*et*np.log(1-xf) - (kap/2+al+et+gg)*np.log(1+xf**2) \
         + lam*(2*gg - 2*xf**2/(1+xf**2)) - 2*kap*chebv(xf, p) - sum(v[k-1]*xf**k/k for k in kk)
    Yf = be*np.log(sf) + gq*np.log(1-sf) + 2*chebv(sf, v)
    val = (-1+al+gg)*np.log(2) + kap*np.sum(p*p/kk) + 0.5*np.sum(v*v/kk) + Xf.max() + Yf.max()
    return val, lam, Xf.max() + Yf.max()


# ---------------------------------------------------------------- objective and descent
def admissible(d):
    n, b, q, g, h = d
    a, L, C, H, A = derived(n, b, q, g, h)
    return b >= 0 and q >= 1 and g >= 1 and h >= 0 and 2*g < n - 1 and 2*b < H

def evaluate(d):
    """quick evaluator: (mu^, gamma^, rho_2^, Sigma_2^, f_low, U_2^)."""
    n, b, q, g, h = d
    fl, _ = f_low(*d)
    v2, _, sig = solve(d, 2)
    gam = fl - float(v2)
    al, gg, gq, et = (b + q)/n, g/n, q/n, h/n
    rho1 = (2 - al - gg)*LOG2 - sig + 2*LOG2
    rho2 = rho1 + (et + gq)*LOG2 + (1 - gg)*LOGS
    mu = rho2*((1 - gg)/gam + n/2) if gam > 0 else float('inf')
    return mu, gam, rho2, float(sig), fl, float(v2)

def descend(n0):
    start = tuple([n0] + [max(1 if i in (1, 2) else 0, round(n0*r)) for i, r in enumerate(START_RATIOS)])
    cache = {}
    def ev(d):
        if d not in cache:
            cache[d] = evaluate(d)
            mu, gam = cache[d][0], cache[d][1]
            print(f"[{time.strftime('%H:%M:%S')}]   n0={n0}: eval {d}: gamma^ {gam:+.6f}, mu^ {mu:.1f}", flush=True)
        return cache[d][0]
    if not admissible(start):
        return {"n0": n0, "start": start, "opt": None, "path": [], "evaluated": {}}
    cur = start; best = ev(cur); path = [(cur, best)]
    while len(cache) < MAXEV:
        cand = []
        for i in range(1, 5):
            for s in (-1, 1):
                d = list(cur); d[i] += s; d = tuple(d)
                if admissible(d): cand.append((ev(d), d))
        if not cand: break
        mu, d = min(cand)
        if mu < best - 1e-3:
            cur, best = d, mu; path.append((cur, best))
        else:
            break
    return {"n0": n0, "start": start, "opt": cur, "mu": best, "values": cache[cur], "path": path,
            "evaluated": {str(k): v for k, v in cache.items()}}

if __name__ == "__main__":
    open(LOGPATH, 'w').close()
    log(f"cat_mu_search.py: N0_LIST {N0_LIST}, WORKERS {WORKERS}, MAXEV {MAXEV}")
    s0 = evaluate((48, 7, 3, 4, 2))[1]
    gate("S1", abs(s0 - 0.001676) < 5e-5, f"start score {s0:+.6f} vs +0.001676")
    log(f"S1 (48, 7, 3, 4, 2): gamma^ {s0:+.6f} (record +0.001676)")
    ctrl = {}
    for d in CONTROLS:
        ctrl[d] = evaluate(d)
        log(f"control {d}: mu^ {ctrl[d][0]:.1f}, gamma^ {ctrl[d][1]:+.6f}, rho_2^ {ctrl[d][2]:.5f}, Sigma_2^ {ctrl[d][3]:.6f}")
    gate("S2", abs(ctrl[CONTROLS[0]][1] - 0.002443) < 5e-5, f"control gamma^ {ctrl[CONTROLS[0]][1]:+.6f} vs +0.002443")
    res = []
    with mp.Pool(WORKERS) as pool:
        for r in pool.imap_unordered(descend, N0_LIST):
            res.append(r)
            if r["opt"] is None:
                log(f"n0={r['n0']}: start {r['start']} not admissible, skipped")
            else:
                v = r["values"]
                log(f"n0={r['n0']}: OPTIMUM {r['opt']}: mu^ {v[0]:.1f}, gamma^ {v[1]:+.6f}, rho_2^ {v[2]:.5f} "
                    f"({len(r['evaluated'])} evaluations; {len(res)}/{len(N0_LIST)} done)")
    log("SUMMARY (ranked by mu^; quick evaluator, not certified)")
    ranked = sorted([r for r in res if r["opt"] is not None], key=lambda r: r["mu"])
    for r in ranked:
        v = r["values"]
        log(f"  {str(r['opt']):22s} mu^ {v[0]:8.1f}   gamma^ {v[1]:+.6f}   rho_2^ {v[2]:.5f}")
    for d in CONTROLS:
        log(f"  control {str(d):14s} mu^ {ctrl[d][0]:8.1f}   gamma^ {ctrl[d][1]:+.6f}")
    total = sum(FAILS.values())
    out = {"script": "cat_mu_search.py", "gates": FAILS, "GATE_FAILS": total,
           "controls": {str(d): v for d, v in ctrl.items()},
           "ranking": [{"design": r["opt"], "values": r["values"]} for r in ranked],
           "runs": res, "seconds": round(time.time() - T0, 1)}
    with open(os.path.join(HERE, "cat_mu_search_out.json"), "w") as fh: json.dump(out, fh, indent=1)
    log(f"gates: {FAILS}")
    print(f"GATE_FAILS = {total}")
    with open(LOGPATH, 'a') as fh: fh.write(f"GATE_FAILS = {total}\n")
    print(f"log saved to {LOGPATH}")
