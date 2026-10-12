"""fz5_check.py -- Fauzan's Hankel determinant for zeta(5) at small sizes (positive control for the beta(4) work).

Run:  python3 fz5_check.py        (no arguments; every switch is set below)
Needs: python-flint (pip install python-flint), mpmath.
Writes fz5_check_log_<stamp>.txt and fz5_check_out_<stamp>.json next to this file.
Expected runtime with 2 workers: about 1-2 minutes for n = 1, 2; about 10-15 minutes with n = 3.

Construction (A. Fauzan, "zeta(5) is irrational", Zenodo, 17 Sep 2026, Section 2): K = 40n, N = 3n, h = 37n,
D_m(t) = prod_{j<=m} (t + j^2), mu_X(t^e) = (-1)^e B_{2e+2} (2e+3)(2e+4)(2e+5)/24,
mu_X(1/(t+j^2)) = j^4 (X - H_j^(5)) - 1/4 + 1/(2j);  G_K(X) = [mu_X(D_N(t)^6 t^(i+j) / D_K(t))]_{i,j<h},
Delta_K = det G_K (polynomial of degree h in X), S_K = (K!)^{2h} 4^{h-1} / ((N!)^{12h} prod_{i<h} ((2i)!)^2), F_K = S_K Delta_K.
Measured: Q = F_K / content(F_K), the primitive integer polynomial (the best possible normalization, at least as good
as the paper's m_{K,M}); log Q(zeta(5)); per n^2 (paper's asymptotic bound: below -139/5 = -27.8); and the bound
b < exp(-log Q(zeta(5))/h): Q excludes zeta(5) = a/b for every such b.
Gates: G1 moments and pole values of mu against the integral with the weight w(y) = (2pi)^4 y^5/12 sum_l l^4 e^{-2 pi l y}
(quadrature, 30 digits, < 1e-25); G2 Q has integer coefficients and degree h, and its leading coefficient is
+-(S_K/content) prod_{N<j<=K} j^4 D_N(-j^2)^5 (paper (2.9)); G3 Q(zeta(5)) > 0 and the same value at two precisions
(relative difference < 1e-30; precision doubled until they agree); G4 (n = 1) Q(zeta(5)) equals the numerical determinant of G_K(zeta(5)) times S_K/content.
"""
import os, sys, math, time, json, datetime
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
from multiprocessing import Pool
import mpmath as mp
from flint import fmpq, fmpq_poly, fmpz_mat

# ---------------- switches ----------------
WORKERS = 2
NS = [1, 2, 3]
# ------------------------------------------

STAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
LOGF = os.path.join(HERE, f"fz5_check_log_{STAMP}.txt"); OUTF = os.path.join(HERE, f"fz5_check_out_{STAMP}.json")
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

def bern(m):
    a, b = mp.bernfrac(m); return fmpq(int(a), int(b))
def mu_t(e): return (-1)**e * bern(2*e+2) * (2*e+3)*(2*e+4)*(2*e+5) / 24
def q2mp(q): return mp.mpf(int(q.p)) / int(q.q)

def entries(n):
    K, N, h = 40*n, 3*n, 37*n
    t = fmpq_poly([0, 1])
    DN = fmpq_poly([1])
    for j in range(1, N+1): DN = DN*(t + j*j)
    Dt = fmpq_poly([1])
    for j in range(N+1, K+1): Dt = Dt*(t + j*j)
    Dtp = Dt.derivative()
    H5 = [fmpq(0)]
    for j in range(1, K+1): H5.append(H5[-1] + fmpq(1, j**5))
    mut = [mu_t(e) for e in range(5*N + 2*h)]
    DN5 = DN**5
    a, b = [], []
    for s in range(2*h - 1):
        q, r = divmod(DN5 * t**s, Dt)
        A = sum((q[e]*mut[e] for e in range(q.degree()+1)), fmpq(0)); B = fmpq(0)
        for j in range(N+1, K+1):
            c = r(-j*j) / Dtp(-j*j)
            B += c * j**4
            A += c * (-j**4 * H5[j] - fmpq(1, 4) + fmpq(1, 2*j))
        a.append(A); b.append(B)
    lead = fmpq(1)
    for j in range(N+1, K+1): lead *= fmpq(j**4) * DN(-j*j)**5
    return K, N, h, a, b, lead

_AI = _BI = None; _H = 0
def _init(ai, bi, h):
    global _AI, _BI, _H
    _AI, _BI, _H = ai, bi, h
def _det(x):
    h = _H
    return x, int(fmpz_mat(h, h, [_AI[i+j] + x*_BI[i+j] for i in range(h) for j in range(h)]).det())

def run(n):
    K, N, h, a, b, lead = entries(n)
    say(f"n = {n}: K = {K}, N = {N}, h = {h}; entries built")
    L = 1
    for q in a + b: L = L*int(q.q)//math.gcd(L, int(q.q))
    ai = [int((q*L).p) for q in a]; bi = [int((q*L).p) for q in b]
    xs = list(range(h+1)); ys = {}
    with Pool(WORKERS, initializer=_init, initargs=(ai, bi, h)) as pool:
        for k, (x, d) in enumerate(pool.imap_unordered(_det, xs)):
            ys[x] = d
            if (k+1) % max(1, (h+1)//4) == 0: say(f"   {k+1}/{h+1} determinants")
    coef = [fmpq(ys[x]) for x in xs]
    for k in range(1, h+1):
        for i in range(h, k-1, -1):
            coef[i] = (coef[i] - coef[i-1]) / (xs[i] - xs[i-k])
    t = fmpq_poly([0, 1])
    P = fmpq_poly([coef[h]])
    for k in range(h-1, -1, -1): P = P*(t - xs[k]) + coef[k]
    P = P / fmpq(L)**h                       # Delta_K(X)
    fac = math.factorial
    SK = fmpq(fac(K)**(2*h) * 4**(h-1), fac(N)**(12*h) * math.prod(fac(2*i)**2 for i in range(1, h)))
    F = P * SK
    g = 0; D = 1
    for k in range(h+1):
        g = math.gcd(g, int(F[k].p)); D = D*int(F[k].q)//math.gcd(D, int(F[k].q))
    content = fmpq(g, D)
    Q = F / content
    okint = Q.degree() == h and all(Q[k].q == 1 for k in range(h+1))
    oklead = abs(Q[h]) == abs(SK / content * lead)
    gate("G2", okint and oklead, f"n={n} integer={okint} leading={oklead}")
    bits = max(abs(int(Q[k].p)).bit_length() for k in range(h+1))
    # the value is far below the coefficients: precision = digits of the largest coefficient + digits of 1/value + margin;
    # start with 1.3 x the coefficient digits and double until two precisions (d and d + 1000) agree to 1e-30
    d = int(1.3*0.302*bits) + 300
    for attempt in range(5):
        vals = []
        for dd in (d, d + 1000):
            mp.mp.dps = dd
            z5 = mp.zeta(5)
            v = mp.mpf(0)
            for k in range(h, -1, -1): v = v*z5 + int(Q[k].p)
            vals.append(v)
        rel = abs(vals[0]-vals[1]) / abs(vals[1])
        if rel < mp.mpf(10)**-30: break
        d *= 2
    gate("G3", vals[1] > 0 and rel < mp.mpf(10)**-30, f"n={n} value {mp.nstr(vals[1], 5)} rel {mp.nstr(rel, 3)}")
    say(f"   value at {d} and {d+1000} digits: relative difference {mp.nstr(rel, 3)}")
    v = vals[1]; logQ = mp.log(v)
    if n == 1:
        G = mp.matrix(h, h)
        for i in range(h):
            for j in range(h): G[i, j] = q2mp(a[i+j]) + mp.zeta(5)*q2mp(b[i+j])
        d2 = mp.det(G) * q2mp(SK/content)
        r4 = abs(d2 - v)/abs(v)
        gate("G4", r4 < mp.mpf(10)**-50, f"rel {mp.nstr(r4, 3)}")
        say(f"   G4 numerical determinant: relative difference {mp.nstr(r4, 3)}")
    lc = math.log(abs(int(content.p))) - math.log(int(content.q))
    ls = math.log(int(SK.p)) - math.log(int(SK.q))
    res = dict(n=n, K=K, N=N, h=h, logQ=float(logQ), logQ_per_n2=float(logQ/n**2), logQ_per_K2=float(logQ/K**2),
               log_content=lc, log_SK=ls, log_max_coeff=bits*math.log(2), log_b_bound=float(-logQ/h))
    say(f"   Q(zeta(5)) = exp({float(logQ):.4f}) > 0;  per n^2: {float(logQ)/n**2:.4f} (paper's asymptotic bound -27.8); per K^2: {float(logQ)/K**2:.5f}")
    say(f"   log content(F_K) = {lc:.2f}, log S_K = {ls:.2f}, log max|coefficient of Q| = {bits*math.log(2):.1f}")
    say(f"   Q excludes zeta(5) = a/b for all b < exp({float(-logQ/h):.4f}) = {float(mp.e**(-logQ/h)):.6g}")
    return res

if __name__ == "__main__":
    say(f"fz5_check.py  workers={WORKERS}  n in {NS}")
    mp.mp.dps = 30
    def w(y):
        return (2*mp.pi)**4 * y**5 / 12 * mp.polylog(-4, mp.e**(-2*mp.pi*y))
    worst = mp.mpf(0)
    for e in range(4):
        I = mp.quad(lambda y: y**(2*e)*w(y), [0, 1, 5, 20, mp.inf]); worst = max(worst, abs(I - q2mp(mu_t(e))))
    for j in (1, 2, 3, 5):
        H = sum(fmpq(1, v**5) for v in range(1, j+1))
        val = j**4*(mp.zeta(5) - q2mp(H)) - mp.mpf(1)/4 + mp.mpf(1)/(2*j)
        I = mp.quad(lambda y: w(y)/(y**2 + j**2), [0, 1, 5, 20, mp.inf]); worst = max(worst, abs(I - val))
    gate("G1", worst < mp.mpf(10)**-25, f"worst {mp.nstr(worst, 3)}")
    say(f"G1 moments e = 0..3 and pole values j = 1, 2, 3, 5 against quadrature: worst difference {mp.nstr(worst, 3)}")
    out = [run(n) for n in NS]
    say("summary  n   K    h   log Q(zeta5)   per n^2    per K^2    log b-bound")
    for r in out:
        say(f"        {r['n']:2d} {r['K']:4d} {r['h']:4d} {r['logQ']:13.3f} {r['logQ_per_n2']:10.3f} {r['logQ_per_K2']:10.5f} {r['log_b_bound']:10.4f}")
    total = sum(FAILS.values())
    json.dump(dict(stamp=STAMP, gates=FAILS, results=out, gate_fails=total), open(OUTF, "w"), indent=1)
    say(f"gates: {FAILS}")
    say(f"GATE_FAILS = {total}")
    say(f"log saved to {LOGF}")
    say(f"json saved to {OUTF}")
