"""b4_fauzan.py -- logged run for Section 6.10 (determinants with poles at the negative squares), 8 October 2026.

Run:  python3 b4_fauzan.py        (no arguments; every switch is set below)
Needs: python-flint, mpmath; full_cat.py, scout_b4n.py, scout_b4l.py, nu_arrays.py in the same folder.
Writes b4_fauzan_log_<stamp>.txt and b4_fauzan_out_<stamp>.json next to this file.
Expected runtime with 2 workers: about 1 minute.

Quantity: Lambda = log|Q(C)|/h^2, Q = Delta/c with c > 0 the rational number that makes Q primitive in Z[X], Delta the
determinant as a polynomial in the constant X, h the size of the matrix, C the constant.
Functionals (Definition fz of the draft): (Z) zeta(s): poles at -(1+qj)^2, mu(t^e) = (2e+s)!|B_{2e+2}|/((s-1)!(2e+2)!),
mu(1/(t+m^2)) = m^(s-1)(X - H_m^(s)) + 1/(2m) - 1/(s-1);  (B) beta(s), s even: poles at -(j+1/2)^2,
mu(t^e) = (2e+s)!(2^(2e+2)-1)|B_{2e+2}|/((s-1)!(2e+1)!(2e+2)), mu(1/(t+(j+1/2)^2)) = 2(2j+1)^(s-1)(-1)^j(X - sum_{l<j}(-1)^l(2l+1)^-s) - 1/(2j+1).
Matrix [mu(D_N(t)^r t^(i+k)/D_K(t))]_{i,k<h}. Arcsine family: OAI (full_cat.py, ratios 48:7:4:4:2), C7 and B7 (scout_b4n.py).

Gates: G1 (B) moments e = 0..3 and pole values j = 0..3 against the integral with w_s(y) = (-1)^(s-1) y^s g^(s-1)(y)/(s-1)!,
g = 1/sinh(pi y), s = 2, 4, 6 (30 digits, < 1e-24); G2 (Z) the same against the integral with 2(-1)^(s-1) y^s f^(s-1)(y)/(s-1)!,
f = 1/(e^(2 pi y) - 1), s = 3, 5, 7; G3 for every item: Q(C) > 0 (Q(C) != 0 for the arcsine family) and two precisions agree
to 1e-30; G4 Fauzan's n = 1 value: log Q(zeta(5)) = -265.129 (fz5_check log _150759) to 1e-3; G5 the arcsine items reproduce the
gaps of b4_route3 (log _110012, Observation nubalance): B7 n = 10, 12: 2.75, 2.94; C7 n = 10, 12: 0.04, 0.10 (to 0.006).
"""
import os, sys, math, time, json, datetime
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE); sys.path.insert(0, HERE)
from multiprocessing import Pool
from fractions import Fraction as Fr
import mpmath as mp
from flint import fmpq, fmpq_poly, fmpz_mat

# ---------------- switches ----------------
WORKERS = 2
# items: (family, K, N, h, r, s, q) for the functionals, ('arc', name, n) for the arcsine family
ITEMS = []
for K in (12, 24):
    for s in (3, 5, 7): ITEMS.append(('Z', K, 0, K, 1, s, 1))
    for s in (2, 4, 6): ITEMS.append(('B', K, 0, K, 1, s, 1))
ITEMS += [('B', 40, 3, 37, 5, 4, 1), ('Z', 40, 3, 37, 6, 5, 1), ('Z', 40, 3, 37, 6, 5, 2)]
for q in (2, 3):
    for s in (3, 5): ITEMS.append(('Z', 24, 0, 24, 1, s, q))
GRID_B4 = [('B', 24, N, h, r, 4, 1) for N in (0, 3, 6, 9) for r in (1, 4, 8) for h in (6, 12, 18) if h <= 24 - N]
GRID_B2 = [('B', 24, N, h, r, 2, 1) for N in (0, 2, 4, 6) for r in (1, 3, 5) for h in (12, 18, 24) if h <= 24 - N]
ARC = [('arc', 'OAI', n) for n in (12, 20, 28, 40)] + [('arc', 'C7', n) for n in (6, 10, 12)] + [('arc', 'B7', n) for n in (6, 10, 12)]
# ------------------------------------------

STAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
LOGF = os.path.join(HERE, f"b4_fauzan_log_{STAMP}.txt"); OUTF = os.path.join(HERE, f"b4_fauzan_out_{STAMP}.json")
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
def q2mp(x): return mp.mpf(int(x.p)) / int(x.q)
def mu_Z(e, s): return fmpq(math.factorial(2*e+s), math.factorial(s-1)*math.factorial(2*e+2)) * abs(bern(2*e+2))
def mu_B(e, s): return fmpq(math.factorial(2*e+s)*(2**(2*e+2)-1), math.factorial(s-1)*math.factorial(2*e+1)*(2*e+2)) * abs(bern(2*e+2))
def const(fam, s, dps):
    mp.mp.dps = dps
    return mp.zeta(s) if fam == 'Z' else (mp.zeta(s, mp.mpf(1)/4) - mp.zeta(s, mp.mpf(3)/4)) / 4**s

def poly_from_matrices(Ai, Bi, h):
    """exact det(Ai + X Bi) (integer matrices) as fmpq_poly, by evaluation at X = 0..h and Newton interpolation"""
    xs = list(range(h+1))
    ys = [fmpq(int(fmpz_mat(h, h, [Ai[i][k] + x*Bi[i][k] for i in range(h) for k in range(h)]).det())) for x in xs]
    coef = ys[:]
    for k in range(1, h+1):
        for i in range(h, k-1, -1): coef[i] = (coef[i]-coef[i-1]) / (xs[i]-xs[i-k])
    t = fmpq_poly([0, 1]); P = fmpq_poly([coef[h]])
    for k in range(h-1, -1, -1): P = P*(t - xs[k]) + coef[k]
    return P
def lam_of(P, cfun, h):
    deg = P.degree()
    g = 0; D = 1
    for k in range(deg+1): g = math.gcd(g, int(P[k].p)); D = D*int(P[k].q)//math.gcd(D, int(P[k].q))
    Q = P * fmpq(D, g)
    bits = max(abs(int(Q[k].p)).bit_length() for k in range(deg+1))
    d = int(1.3*0.302*bits) + 300
    for _ in range(5):
        vals = []
        for dd in (d, d+1000):
            C = cfun(dd); v = mp.mpf(0)
            for k in range(deg, -1, -1): v = v*C + int(Q[k].p)
            vals.append(v)
        rel = abs(vals[0]-vals[1]) / abs(vals[1])
        if rel < mp.mpf(10)**-30: break
        d *= 2
    return float(mp.log(abs(vals[1]))), vals[1] > 0, float(rel) if rel > 0 else 0.0, deg

def fz_item(fam, K, N, h, r, s, q):
    t = fmpq_poly([0, 1])
    if fam == 'Z':
        a2 = [fmpq((1+q*j)**2) for j in range(K)]
        mmax = 1 + q*(K-1); H = [fmpq(0)]
        for v in range(1, mmax+1): H.append(H[-1] + fmpq(1, v**s))
        def pole(j):
            m = 1 + q*j; return fmpq(m**(s-1)), -fmpq(m**(s-1))*H[m] + fmpq(1, 2*m) - fmpq(1, s-1)
        mu = mu_Z
    else:
        a2 = [fmpq((2*j+1)**2, 4) for j in range(K)]
        Bp = [fmpq(0)]
        for l in range(K): Bp.append(Bp[-1] + fmpq((-1)**l, (2*l+1)**s))
        def pole(j):
            k = 2*(2*j+1)**(s-1)*(-1)**j; return fmpq(k), -k*Bp[j] - fmpq(1, 2*j+1)
        mu = mu_B
    DN = fmpq_poly([1])
    for j in range(N): DN = DN*(t + a2[j])
    Dt = fmpq_poly([1])
    for j in range(N, K): Dt = Dt*(t + a2[j])
    Dtp = Dt.derivative(); num0 = DN**(r-1)
    mus = [mu(e, s) for e in range((r-1)*N + 2*h)]
    PV = [pole(j) for j in range(K)]
    A, B = [], []
    for e in range(2*h-1):
        qq, rem = divmod(num0 * t**e, Dt)
        a = sum((qq[k]*mus[k] for k in range(qq.degree()+1)), fmpq(0)); b = fmpq(0)
        for j in range(N, K):
            c = rem(-a2[j]) / Dtp(-a2[j]); bx, ax = PV[j]
            b += c*bx; a += c*ax
        A.append(a); B.append(b)
    L = 1
    for x in A + B: L = L*int(x.q)//math.gcd(L, int(x.q))
    Ai = [[int((A[i+k]*L).p) for k in range(h)] for i in range(h)]; Bi = [[int((B[i+k]*L).p) for k in range(h)] for i in range(h)]
    P = poly_from_matrices(Ai, Bi, h)
    lq, pos, rel, deg = lam_of(P, lambda dd: const(fam, s, dd), h)
    return lq, pos, rel, deg

def arc_item(name, n):
    from nu_arrays import c
    if name == 'OAI':
        import full_cat as FC
        A, B = FC.build(n, round(7*n/48), round(4*n/48), round(4*n/48), round(2*n/48)); key = 'G'
    else:
        import scout_b4n as Nm
        Mnu, MA, ML = Nm.arrays(2*n+6)
        if name == 'C7':
            E = [[MA(2*r, 2*k) for k in range(n)] for r in range(n)]
            A = [[e.get('1', Fr(0)) for e in row] for row in E]; B = [[e.get('G', Fr(0)) for e in row] for row in E]; key = 'G'
        else:
            J = [2*k for k in range(n)]; A = [[None]*n for _ in range(n)]; B = [[None]*n for _ in range(n)]; key = 'b4'
            for r in range(n):
                col = {j: Nm.addd({}, Mnu(2*r, j)) for j in J}
                aL = {}
                for j in sorted(J, reverse=True): aL[j] = -col[j].get('X3', Fr(0)) - sum((aL[i2]*c((i2-j)//2) for i2 in aL if i2 > j), Fr(0))
                for j in J:
                    for i2, av in aL.items(): col[j] = Nm.addd(col[j], ML(i2, j), av)
                aA = {}
                for j in sorted(J, reverse=True): aA[j] = (-col[j].get('G', Fr(0)) - sum((aA[i2]*4*c((i2-j)//2) for i2 in aA if i2 > j), Fr(0)))/4
                for j in J:
                    for i2, av in aA.items(): col[j] = Nm.addd(col[j], MA(i2, j), av)
                for k in range(n):
                    assert all(kk in ('1', 'b4') for kk in col[2*k]); A[r][k] = col[2*k].get('1', Fr(0)); B[r][k] = col[2*k].get('b4', Fr(0))
    L = 1
    for row in A + B:
        for v in row: L = L*v.denominator//math.gcd(L, v.denominator)
    Ai = [[int(v*L) for v in row] for row in A]; Bi = [[int(v*L) for v in row] for row in B]
    P = poly_from_matrices(Ai, Bi, n)
    cf = (lambda dd: (mp.mp.__setattr__('dps', dd), +mp.catalan)[1]) if key == 'G' else (lambda dd: const('B', 4, dd))
    lq, pos, rel, deg = lam_of(P, cf, n)
    return lq, True, rel, deg

def work(it):
    t = time.time()
    if it[0] == 'arc':
        lq, pos, rel, deg = arc_item(it[1], it[2]); h = it[2]
    else:
        lq, pos, rel, deg = fz_item(*it); h = it[3]
    return dict(item=list(it), logQ=lq, Lambda=lq/h**2, perK2=(lq/it[1]**2 if it[0] != 'arc' else None),
                positive=bool(pos), rel=rel, deg=deg, secs=time.time()-t)

def name(it):
    if it[0] == 'arc': return f"{it[1]} n={it[2]}"
    fam, K, N, h, r, s, q = it
    return f"{'zeta' if fam == 'Z' else 'beta'}({s}) K={K} N={N} h={h} r={r}" + (f" q={q}" if q != 1 else "")

if __name__ == "__main__":
    say(f"b4_fauzan.py  workers={WORKERS}  items={len(ITEMS)} + grids {len(GRID_B4)}, {len(GRID_B2)} + arcsine {len(ARC)}")
    mp.mp.dps = 30
    wB = lambda s: (lambda y: (-1)**(s-1)*y**s*mp.diff(lambda u: 1/mp.sinh(mp.pi*u), y, s-1)/math.factorial(s-1))
    worst = mp.mpf(0)
    for s in (2, 4, 6):
        w = wB(s); b = const('B', s, 30)
        for e in range(4):
            I = mp.quad(lambda y: y**(2*e)*w(y), [0, 0.5, 2, 8, 30, mp.inf]); worst = max(worst, abs(I - q2mp(mu_B(e, s))))
        for j in range(4):
            a = mp.mpf(2*j+1)/2; Bj = sum(mp.mpf((-1)**l)/(2*l+1)**s for l in range(j))
            val = 2*(2*j+1)**(s-1)*(-1)**j*(b - Bj) - mp.mpf(1)/(2*j+1)
            I = mp.quad(lambda y: w(y)/(y**2 + a**2), [0, 0.5, 2, 8, 30, mp.inf]); worst = max(worst, abs(I - val))
    gate("G1", worst < mp.mpf(10)**-24, f"worst {mp.nstr(worst, 3)}")
    say(f"G1 functional (B), s = 2, 4, 6, against its weight: worst difference {mp.nstr(worst, 3)}")
    worst = mp.mpf(0)
    for s in (3, 5, 7):
        w = (lambda s: lambda y: 2*(-1)**(s-1)*y**s*mp.diff(lambda u: 1/(mp.e**(2*mp.pi*u)-1), y, s-1)/math.factorial(s-1))(s)
        z = mp.zeta(s)
        for e in range(4):
            I = mp.quad(lambda y: y**(2*e)*w(y), [0, 0.5, 2, 8, 30, mp.inf]); worst = max(worst, abs(I - q2mp(mu_Z(e, s))))
        for m in (1, 2, 3):
            Hm = sum(mp.mpf(1)/v**s for v in range(1, m+1))
            val = m**(s-1)*(z - Hm) + mp.mpf(1)/(2*m) - mp.mpf(1)/(s-1)
            I = mp.quad(lambda y: w(y)/(y**2 + m**2), [0, 0.5, 2, 8, 30, mp.inf]); worst = max(worst, abs(I - val))
    gate("G2", worst < mp.mpf(10)**-24, f"worst {mp.nstr(worst, 3)}")
    say(f"G2 functional (Z), s = 3, 5, 7, against its weight: worst difference {mp.nstr(worst, 3)}")
    allitems = ITEMS + GRID_B4 + GRID_B2 + ARC
    res = {}
    with Pool(WORKERS) as pool:
        for k, o in enumerate(pool.imap_unordered(work, allitems)):
            it = tuple(o['item']); res[it] = o
            gate("G3", (o['positive'] or it[0] == 'arc') and o['rel'] < 1e-30, f"{name(it)} positive={o['positive']} rel={o['rel']}")
            if (k+1) % 10 == 0: say(f"   {k+1}/{len(allitems)} items")
    def R(it): return res[tuple(it)]
    say("Observation fzweights (N = 0, r = 1, h = K): Lambda")
    for K in (12, 24):
        say(f"   K = {K}: " + "  ".join(f"zeta({s}) {R(('Z', K, 0, K, 1, s, 1))['Lambda']:+.3f}" for s in (3, 5, 7)) + "  " +
            "  ".join(f"beta({s}) {R(('B', K, 0, K, 1, s, 1))['Lambda']:+.3f}" for s in (2, 4, 6)))
    o = R(('B', 40, 3, 37, 5, 4, 1))
    say(f"   beta(4), Fauzan's parameters n = 1 (K 40, N 3, h 37, r 5): Lambda {o['Lambda']:+.3f}, per K^2 {o['perK2']:+.3f}")
    g4 = [R(it) for it in GRID_B4]
    say(f"   beta(4) grid, {len(g4)} designs at K = 24: smallest log Q/K^2 {min(x['perK2'] for x in g4):+.3f}, smallest Lambda {min(x['Lambda'] for x in g4):+.3f}")
    say("Observation fzstep (N = 0, r = 1, h = K = 24): Lambda for pole step q = 1, 2, 3")
    for s in (3, 5):
        say(f"   zeta({s}): " + "  ".join(f"q={q} {R(('Z', 24, 0, 24, 1, s, q))['Lambda']:+.3f}" for q in (1, 2, 3)))
    a1 = R(('Z', 40, 3, 37, 6, 5, 1)); a2 = R(('Z', 40, 3, 37, 6, 5, 2))
    say(f"   zeta(5), Fauzan's parameters n = 1: log Q/K^2 q=1 {a1['perK2']:+.4f}, q=2 {a2['perK2']:+.4f}")
    gate("G4", abs(a1['logQ'] + 265.129) < 1e-3, f"logQ {a1['logQ']}")
    say("Observation fzscale: Lambda on the common scale")
    for nm, ns in (('OAI', (12, 20, 28, 40)), ('C7', (6, 10, 12)), ('B7', (6, 10, 12))):
        say(f"   {nm}: " + "  ".join(f"{n}: {R(('arc', nm, n))['Lambda']:+.4f}" for n in ns))
    say(f"   beta(2) (B): 12: {R(('B', 12, 0, 12, 1, 2, 1))['Lambda']:+.4f}  24: {R(('B', 24, 0, 24, 1, 2, 1))['Lambda']:+.4f}")
    say(f"   beta(4) (B): 12: {R(('B', 12, 0, 12, 1, 4, 1))['Lambda']:+.4f}  24: {R(('B', 24, 0, 24, 1, 4, 1))['Lambda']:+.4f}  37: {o['Lambda']:+.4f}")
    g2 = [R(it) for it in GRID_B2]
    say(f"   beta(2) grid, {len(g2)} designs at K = 24: smallest Lambda {min(x['Lambda'] for x in g2):+.4f}, all positive: {all(x['Lambda'] > 0 for x in g2)}")
    ref = {('B7', 10): 2.75, ('B7', 12): 2.94, ('C7', 10): 0.04, ('C7', 12): 0.10}
    ok5 = all(abs(R(('arc', k[0], k[1]))['Lambda'] - v) < 0.006 for k, v in ref.items())
    gate("G5", ok5, str({k: R(('arc', k[0], k[1]))['Lambda'] for k in ref}))
    total = sum(FAILS.values())
    json.dump(dict(stamp=STAMP, gates=FAILS, results=[res[tuple(it)] for it in allitems], gate_fails=total), open(OUTF, "w"), indent=1)
    say(f"gates: {FAILS}")
    say(f"GATE_FAILS = {total}")
    say(f"log saved to {LOGF}")
    say(f"json saved to {OUTF}")
