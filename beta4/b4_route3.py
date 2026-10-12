"""b4_route3.py -- logged runs for the beta(4) constructions of the route-3 work (8 October 2026).

Run:  python3 b4_route3.py        (no arguments; every switch is set below)
Needs, in the same folder: nu_arrays.py, scout_b4l.py, scout_b4n.py, full_cat.py, rz_scout.py, rz_hankel.py, b4_balance.py.
Writes b4_route3_log_<stamp>.txt and b4_route3_out_<stamp>.json next to this file.
Expected runtime with 2 workers: about 4-6 minutes (the beta(4) overlap determinant at n = 32 dominates).

Quantities (as in b4_balance.py): lambda = log|Delta|/n^2 with the true constant; varpi = -log den(Delta with the
constant(s) replaced by 1/3)/n^2; gap = lambda - varpi.

R1  positive controls (gates):
    C1  exact arrays (new measure nu, arcsine kernel 1/(1-ts), arcsine kernel (-log s)/(1-ts)) against quadrature,
        6 index pairs each, 30 digits: |difference| < 1e-25.
    C2  the preprint's determinant at its own design (48,7,4,4,2): lambda - log2/2 = -2.2356 (cat_extend.py record), to 1e-4.
    C3  Rivoal-Zudilin moments, exact against direct summation (eta = (3,1,1,1,1,1,1), N = 2, j = 0, 2): < 1e-15.
R2  beta(4) overlap design (#7: measure nu, even indices, X3 and G cancelled) and Catalan with the same rows (#7c):
    n = 4, 6, 8, 10, 12, 16, 20, 24, 28, 32 with h = q = 0; and h = q = n/4 for n = 8, 12, 16.
    Gates: no constant other than beta(4) left in any entry (C4); two precisions agree to 1e-30 (C5); determinants nonzero (C6).
R3  the preprint's determinant for G at ratios 48:7:4:4:2, n = 12, 20, 28, 48.
R4  denominator structure at n = 16: per-prime valuations (p <= 13) with row bound and largest entry power, and the
    split of log den / n^2 into bins of p/n, for #7, #7c and the preprint's determinant.
R5  Rivoal-Zudilin Hankel determinants: B = log|H_m| + log den(H_m with beta(2), beta(4) -> 1/3), eta = (3,1^5,1) N = 12,
    m = 1, 2, 3, 4, 6, 8, 10; eta = (5,2^5,2) N = 6, m = 1, 2; and the maximal m for N = 3, 4, 6, 8, 10, 12, 14.
"""
import os, sys, math, time, json, datetime
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE); sys.path.insert(0, HERE)
os.environ.setdefault("OMP_NUM_THREADS", "2")
from fractions import Fraction as Fr
from math import comb
from multiprocessing import Pool
import mpmath as mp

# ---------------- switches ----------------
WORKERS = 2
NS_R2 = [4, 6, 8, 10, 12, 16, 20, 24, 28, 32]
NS_R2_HQ = [8, 12, 16]
NS_R3 = [12, 20, 28, 48]
N_R4 = 16
R5_FIXED = [((3, 1, 1, 1, 1, 1, 1), 12, [1, 2, 3, 4, 6, 8, 10]), ((5, 2, 2, 2, 2, 2, 2), 6, [1, 2])]
R5_MAXM_NS = [3, 4, 6, 8, 10, 12, 14]
TOL_C1, TOL_C2, TOL_C3, TOL_PREC = mp.mpf('1e-25'), 1e-4, mp.mpf('1e-15'), mp.mpf('1e-30')
# ------------------------------------------

STAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
LOGF = os.path.join(HERE, f"b4_route3_log_{STAMP}.txt"); OUTF = os.path.join(HERE, f"b4_route3_out_{STAMP}.json")
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

def B7_matrix(n, h=0, q=0):
    """#7 overlap design: rows t^{2r}(1-t^2)^h, columns s^{2k}(1-s^2)^q; returns (Frat, Fb4, residues)."""
    import scout_b4n as N
    from nu_arrays import c
    top = 2*(n-1+max(h, q)) + 4
    Mnu, MA, ML = N.arrays(top)
    cols = [{2*(k+l): Fr((-1)**l*comb(q, l)) for l in range(q+1)} for k in range(n)]
    J = sorted(set(j for cc in cols for j in cc))
    Frat = [[Fr(0)]*n for _ in range(n)]; Fb4 = [[Fr(0)]*n for _ in range(n)]; resid = 0
    for r in range(n):
        a = {2*(r+l): Fr((-1)**l*comb(h, l)) for l in range(h+1)}
        col = {j: {} for j in J}
        for j in J:
            for i, av in a.items(): col[j] = N.addd(col[j], Mnu(i, j), av)
        aL = {}
        for j in sorted(J, reverse=True): aL[j] = -col[j].get('X3', Fr(0)) - sum((aL[i2]*c((i2-j)//2) for i2 in aL if i2 > j), Fr(0))
        for j in J:
            for i2, av in aL.items(): col[j] = N.addd(col[j], ML(i2, j), av)
        aA = {}
        for j in sorted(J, reverse=True): aA[j] = (-col[j].get('G', Fr(0)) - sum((aA[i2]*4*c((i2-j)//2) for i2 in aA if i2 > j), Fr(0)))/4
        for j in J:
            for i2, av in aA.items(): col[j] = N.addd(col[j], MA(i2, j), av)
        for k, ck in enumerate(cols):
            e = {}
            for j, cv in ck.items(): e = N.addd(e, col[j], cv)
            if any(kk not in ('1', 'b4') for kk in e): resid += 1
            Frat[r][k] = e.get('1', Fr(0)); Fb4[r][k] = e.get('b4', Fr(0))
    return Frat, Fb4, resid
def C7_matrix(n):
    import scout_b4n as N
    Mnu, MA, ML = N.arrays(2*n+6)
    E = [[MA(2*r, 2*k) for k in range(n)] for r in range(n)]
    resid = sum(1 for row in E for e in row if any(kk not in ('1', 'G') for kk in e))
    return [[e.get('1', Fr(0)) for e in row] for row in E], [[e.get('G', Fr(0)) for e in row] for row in E], resid
def PRE_matrix(n):
    import full_cat as FC
    return FC.build(n, round(7*n/48), round(4*n/48), round(4*n/48), round(2*n/48))
def measure(Frat, Fc, const):
    from scout_b4l import det_exact
    n = len(Frat)
    Dh = det_exact([[Frat[r][k] + Fr(1, 3)*Fc[r][k] for k in range(n)] for r in range(n)])
    outs = []
    for dps in (60 + 6*n*n, 120 + 6*n*n):
        mp.mp.dps = dps
        C = (mp.zeta(4, mp.mpf(1)/4) - mp.zeta(4, mp.mpf(3)/4))/256 if const == 'b4' else +mp.catalan
        Mm = mp.matrix(n, n)
        for r in range(n):
            for k in range(n):
                Mm[r, k] = mp.mpf(Frat[r][k].numerator)/Frat[r][k].denominator + C*mp.mpf(Fc[r][k].numerator)/Fc[r][k].denominator
        outs.append(mp.det(Mm))
    rel = abs(outs[0]-outs[1])/abs(outs[1]) if outs[1] != 0 else mp.inf
    lam = float(mp.log(abs(outs[1]))/n**2); vp = -math.log(Dh.denominator)/n**2 if Dh != 0 else None
    return dict(n=n, lam=lam, varpi=vp, gap=lam-vp, rel=float(rel), nonzero=bool(outs[1] != 0 and Dh != 0), den=Dh.denominator)
def task(t):
    kind, n, h, q = t; t0 = time.time()
    if kind == 'B7': Fa, Fb, res = B7_matrix(n, h, q); R = measure(Fa, Fb, 'b4')
    elif kind == 'C7': Fa, Fb, res = C7_matrix(n); R = measure(Fa, Fb, 'G')
    else: Fa, Fb = PRE_matrix(n); res = 0; R = measure(Fa, Fb, 'G')
    R.update(kind=kind, h=h, q=q, resid=res, secs=time.time()-t0); R['den'] = str(R['den'])
    return R

def primes(N_):
    s = [1]*(N_+1); s[0] = s[1] = 0
    for i in range(2, int(N_**0.5)+1):
        if s[i]:
            for j in range(i*i, N_+1, i): s[j] = 0
    return [i for i in range(N_+1) if s[i]]
def vpow(x, p):
    v = 0
    while x % p == 0: x //= p; v += 1
    return v

if __name__ == "__main__":
    say(f"b4_route3.py  workers={WORKERS}")
    # ---------------- R1 controls
    import scout_b4n as N
    mp.mp.dps = 30
    pi = mp.pi; G = +mp.catalan; X3 = 8*mp.im(mp.polylog(3, mp.mpc(1, 1)/2)) - pi**3/16 - pi/4*mp.log(2)**2
    B4 = (mp.zeta(4, mp.mpf(1)/4) - mp.zeta(4, mp.mpf(3)/4))/256
    val = {'1': 1, 'G': G, 'X3': X3, 'b4': B4}
    Mnu, MA, ML = N.arrays(30)
    def lm(th): return -mp.log(2)-2*mp.log(mp.sin(pi/4-th/2))
    def lp(th): return -mp.log(1+mp.sin(th))
    def ker(j, t, Lg, p):
        if abs(t) < mp.mpf('0.3'): return mp.nsum(lambda u: t**u/(u+j+1)**p, [0, mp.inf])
        head = Lg if p == 1 else mp.polylog(2, t)
        return (head - sum(t**k/mp.mpf(k)**p for k in range(1, j+1)))/t**(j+1)
    def quadv(w, i, j, p):
        f = lambda th: w(th)*(mp.sin(th)**i*ker(j, mp.sin(th), lm(th), p) + (-mp.sin(th))**i*ker(j, -mp.sin(th), lp(th), p))
        return mp.re(mp.quad(f, [0, pi/4, pi/2]))
    wnu = lambda th: mp.sin(th)*(th*th - pi*th); wA = lambda th: mp.sin(th)
    worst = mp.mpf(0)
    for (i, j) in [(0, 0), (0, 2), (2, 6), (6, 2), (8, 8), (2, 10)]:
        for nm, F_, w, p in (('nu', Mnu, wnu, 1), ('A', MA, wA, 1), ('L', ML, wA, 2)):
            ex = sum(mp.mpf(v.numerator)/v.denominator*val[k] for k, v in F_(i, j).items())
            worst = max(worst, abs(ex - quadv(w, i, j, p)))
    gate("C1", worst < TOL_C1, f"worst {mp.nstr(worst, 3)}")
    say(f"C1 exact arrays vs quadrature, 18 cases at 30 digits: worst {mp.nstr(worst, 3)}")
    import rz_hankel as RZ
    mp.mp.dps = 40
    eta = (3, 1, 1, 1, 1, 1, 1); Nn = 2; qq = 6; h0 = 3*Nn+1
    hh = [h0]+[mp.mpf(2*eta[j]*Nn+1)/2 for j in range(1, qq)]+[eta[qq]*Nn+1]
    g = mp.mpf(1)
    for j in range(1, qq): g *= mp.factorial((eta[0]-2*eta[j])*Nn)
    g /= mp.factorial(eta[qq]*Nn)**2
    wr = mp.mpf(0)
    for jj in (0, 2):
        c0, c2, c4, _ = RZ.moment(eta, Nn, jj)
        ex = mp.mpf(c0.numerator)/c0.denominator + mp.mpf(c2.numerator)/c2.denominator*G + mp.mpf(c4.numerator)/c4.denominator*B4
        def term(t):
            v = (2*t+h0)*mp.gamma(t+h0)/mp.gamma(t+1)*(t+mp.mpf(h0)/2)**(2*jj)
            for j in range(1, qq+1): v *= mp.gamma(t+hh[j])/mp.gamma(t+1+h0-hh[j])
            return g*v*(-1)**int(t)
        wr = max(wr, abs(ex - mp.nsum(term, [0, mp.inf])))
    gate("C3", wr < TOL_C3, f"worst {mp.nstr(wr, 3)}")
    say(f"C3 Rivoal-Zudilin moments exact vs direct series: worst {mp.nstr(wr, 3)}")
    # ---------------- R2, R3 (parallel)
    tasks = [('B7', n, 0, 0) for n in NS_R2] + [('B7', n, n//4, n//4) for n in NS_R2_HQ] + [('C7', n, 0, 0) for n in NS_R2] \
            + [('PRE', n, 0, 0) for n in NS_R3]
    tasks.sort(key=lambda t: -t[1])
    results = []
    with Pool(WORKERS) as pool:
        for R in pool.imap_unordered(task, tasks):
            results.append(R)
            if R['kind'] == 'B7': gate("C4", R['resid'] == 0, f"n={R['n']} residues {R['resid']}")
            gate("C5", R['rel'] < TOL_PREC, f"{R['kind']} n={R['n']} rel {R['rel']:.1e}")
            gate("C6", R['nonzero'], f"{R['kind']} n={R['n']}")
            extra = f" lambda-log2/2 = {R['lam']-math.log(2)/2:.5f}" if R['kind'] == 'PRE' else ""
            say(f"{R['kind']:3s} n={R['n']:2d} h={R['h']} q={R['q']}: lambda {R['lam']:.4f} varpi {R['varpi']:.4f} gap {R['gap']:.4f}"
                f"{extra}  prec {R['rel']:.1e}  residues {R['resid']}  ({R['secs']:.1f}s)")
    pre48 = [R for R in results if R['kind'] == 'PRE' and R['n'] == 48][0]
    d = abs((pre48['lam'] - math.log(2)/2) - (-2.2356))
    gate("C2", d < TOL_C2, f"diff {d:.2e}")
    say(f"C2 preprint design (48,7,4,4,2): lambda - log2/2 = {pre48['lam']-math.log(2)/2:.5f} (record -2.2356), diff {d:.1e}")
    # ---------------- R4 denominator structure at n = N_R4
    from scout_b4l import det_exact
    n = N_R4; R4 = {}
    for nm, (Fa, Fb) in (('B7', B7_matrix(n)[:2]), ('C7', C7_matrix(n)[:2]), ('PRE', PRE_matrix(n))):
        M = [[Fa[r][k] + Fr(1, 3)*Fb[r][k] for k in range(n)] for r in range(n)]
        D = det_exact(M).denominator
        small = {}
        for p in (3, 5, 7, 11, 13):
            small[p] = (vpow(D, p), sum(max(vpow(x.denominator, p) for x in row) for row in M), max(vpow(x.denominator, p) for row in M for x in row))
        bins = {}
        rest = D
        for p in primes(8*n+20):
            v = vpow(rest, p)
            if v:
                rest //= p**v; b = min(int(4*p/n), 15); bins[b] = bins.get(b, 0.0) + v*math.log(p)/n**2
        gate("C7", rest == 1, f"{nm} unfactored part")
        R4[nm] = dict(small={str(k): v for k, v in small.items()}, bins={str(b/4): x for b, x in sorted(bins.items())}, total=math.log(D)/n**2)
        say(f"R4 {nm} n={n}: total {math.log(D)/n**2:.3f}; bins p/n: " + " ".join(f"[{b/4:.2f}):{x:.2f}" for b, x in sorted(bins.items())))
        say(f"R4 {nm} n={n}: (det, row bound, largest entry) " + " ".join(f"p={p}:{v}" for p, v in small.items()))
    # ---------------- R5 Rivoal-Zudilin Hankel
    R5 = []
    def hankel(eta, Nn, m):
        mus = [RZ.moment(eta, Nn, j) for j in range(2*m-1)]
        mp.mp.dps = 80 + 3*m*m + 40*Nn
        b2 = +mp.catalan; b4 = (mp.zeta(4, mp.mpf(1)/4)-mp.zeta(4, mp.mpf(3)/4))/256
        fv = lambda x: mp.mpf(x.numerator)/x.denominator
        H = mp.matrix(m, m)
        for r in range(m):
            for k in range(m): c0, c2, c4, _ = mus[r+k]; H[r, k] = fv(c0)+fv(c2)*b2+fv(c4)*b4
        Ht = mp.det(H)
        Dh = det_exact([[mus[r+k][0]+Fr(1, 3)*(mus[r+k][1]+mus[r+k][2]) for k in range(m)] for r in range(m)])
        gate("C6", Ht != 0 and Dh != 0, f"Hankel eta={eta} N={Nn} m={m}")
        return float(mp.log(abs(Ht))), math.log(Dh.denominator)
    for eta, Nn, ms in R5_FIXED:
        for m in ms:
            a, b = hankel(eta, Nn, m); R5.append(dict(eta=eta, N=Nn, m=m, logH=a, logden=b, B=a+b))
            say(f"R5 eta={eta} N={Nn} m={m}: log|H| {a:.1f}  log den {b:.1f}  B {a+b:+.1f}")
    for Nn in R5_MAXM_NS:
        m = (3*Nn+6)//4; a, b = hankel((3, 1, 1, 1, 1, 1, 1), Nn, m)
        R5.append(dict(eta=(3, 1, 1, 1, 1, 1, 1), N=Nn, m=m, logH=a, logden=b, B=a+b, gap_per_m2=(a+b)/m**2))
        say(f"R5 eta=(3,1^5,1) N={Nn} maximal m={m}: B {a+b:+.1f}, per m^2 {(a+b)/m**2:+.3f}")
    total = sum(FAILS.values())
    json.dump(dict(stamp=STAMP, gates=FAILS, R2R3=results, R4=R4, R5=R5, gate_fails=total), open(OUTF, "w"), indent=1, default=str)
    say(f"gates: {FAILS}")
    say(f"GATE_FAILS = {total}")
    say(f"log saved to {LOGF}")
    say(f"json saved to {OUTF}")
