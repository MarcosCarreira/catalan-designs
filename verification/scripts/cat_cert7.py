"""
cat_cert7.py -- independent rebuild of the Section 7 certificate of
"Catalan's constant is irrational" (OpenAI, 24 Sep 2026), pp. 34-41
(Proposition 7.1: limsup (log|Delta_N|/n^2 - log2/2) <= -2.290939875).

INPUT DATA: transcribed by hand from the paper (tables p.35, p.37-40);
all integer coefficients x 1e-8; nonreal base pairs (a,b) mean
r_z = (a - i b)/2 * 1e-8 and r_conj(z) = conj(r_z).

WHAT IT CHECKS (pre-registered; expected results in brackets)
  G1  exact derivative numerators A = Q*X', Q*Y' over Q(i): every polynomial
      division exact, imaginary parts cancel, degrees 36, 24, 13, 9    [0 fails]
  G2  A/Q equals the numerical derivative of X, Y (mpmath, 40 dps) at
      5 points per function, relative error < 1e-25                     [0 fails]
  G3  Descartes sign variations on the division intervals (eq. 94):
      k=2 X: 9,9 | k=2 Y: 5,2,2,6 | k=1 X: 1,1,2,1 | k=1 Y: 5          [0 fails]
  G4  root brackets (m, m+2)/1e10 (p.38): exact sign change of A at both
      ends for all 43, and bracket count per interval = variation count [0 fails]
  G5  Table 2: X, Y evaluated at the 50 listed points at 40 dps (ordinary
      logarithms, independent of the paper's rational-log recipe):
      value <= listed upper bound                                       [0 fails]
      soft: listed - value < 1.1e-12 (bounds "rounded upward to twelve
      decimals" after an error allowance below 1e-13)                   [report]
  G6  norms k*||p||^2 + ||v||^2/2: value <= .778415976284 (k=2),
      .931985203901 (k=1); soft gap < 1.1e-12                           [0 fails]
  G7  cutoffs (98): every Table 2 value and every computed value
      <= -.98399/-1.60890 (k=2), -1.3244/-1.4280 (k=1); norms <= .77844/.9321;
      final sums with log2 > .693146 give exactly -2.290939875 and
      -2.296789875                                                       [0 fails]
  G8  positive control of the global-sup claim: dense grid (40000 points
      per domain, float arithmetic) max of X, Y <= (largest Table 2 value)
      + 0.000024 (the paper's mean-value allowance, Sec. 7.4)            [0 fails]
  G9  derivative bound of Sec. 7.4: |X'|, |Y'| < 12852 at every bracket
      left endpoint (paper: < 12852 throughout each bracket)             [0 fails]

RUN: python3 cat_cert7.py     (no arguments; ~2 s; single process; needs mpmath)
OUTPUT: progress lines on screen AND in <name>_log_<YYYYMMDD_HHMMSS>.txt; cat_cert7_out.json;
        final GATE_FAILS line. All files are written next to the script.
"""
import sys, json, time, cmath, math
from fractions import Fraction as Fr
try:
    sys.set_int_max_str_digits(0)
except AttributeError:
    pass
import mpmath as mpm

T0 = time.time()
import os as _os
HERE = _os.path.dirname(_os.path.abspath(__file__))
STAMP = time.strftime('%Y%m%d_%H%M%S')
LOGPATH = _os.path.join(HERE, 'cat_cert7_log_' + STAMP + '.txt')
LOGF = open(LOGPATH, 'w')
def log(s):
    line = f"[{time.time()-T0:6.1f}s] {s}"
    print(line, flush=True)
    LOGF.write(line + '\n'); LOGF.flush()
FAILS = {}
SOFT = []
def gate(name, ok, info=""):
    FAILS[name] = FAILS.get(name, 0) + (0 if ok else 1)
    if not ok:
        log(f"FAIL {name} {info}")

# ------------------------------------------------------------ data (paper)
ALPHA, BETA, GAMMA, ETA = Fr(11, 48), Fr(7, 48), Fr(4, 48), Fr(2, 48)
E8 = Fr(1, 10**8)
FIN = {  # finite coefficients l_1..l_d
 (2, 'p'): [45559127, -50750856, -6578767, 13970217, 4786184,
            -4292433, -576704, 1311615, 671564, -346453],
 (2, 'v'): [-23910158, 21152432, -2885110, -11558199, 6485289,
            1456821, -1912176, -2263524, 2742210, -1162454],
 (1, 'p'): [11913521, -79993701, -21956443, 37903579,
            10744022, -2073193, 570246, -24939103],
 (1, 'v'): [-89913025, 52874280, 45168341, -30708629,
            -19269841, 33922112, 9819141, -16992389]}
LAM = {2: Fr(0), 1: Fr('2.47405979')}
# exponential terms, k=2 only: (re z, im z, a, b); b None for real base
EXP = {
 (2, 'p'): [('.85', '0', -9338452, None), ('.94', '0', -2141509, None),
            ('0', '.7', -66277922, -31907569), ('0', '.85', -1231651, 6002645),
            ('.092', '.92', -3225918, 8928234), ('-.092', '.92', -2105536, -9091287)],
 (2, 'v'): [('-.8', '0', 15199211, None), ('-.96', '0', 4451662, None),
            ('.88', '0', 2545398, None), ('.95', '0', -4932634, None),
            ('.984', '0', 11618157, None),
            ('0', '.78', 27238714, -38447936), ('0', '.9', -31341084, -30188786),
            ('0', '.955', -6693542, 11912254), ('0', '.984', 2055213, -21715849)],
 (1, 'p'): [], (1, 'v'): []}
DIV = {(2, 'X'): ['-1', '0', '1'], (2, 'Y'): ['0', '.25', '.5', '.75', '1'],
       (1, 'X'): ['-1', '-.5', '0', '.5', '1'], (1, 'Y'): ['0', '1']}
VAR = {(2, 'X'): [9, 9], (2, 'Y'): [5, 2, 2, 6], (1, 'X'): [1, 1, 2, 1], (1, 'Y'): [5]}
DEG = {(2, 'X'): 36, (2, 'Y'): 24, (1, 'X'): 13, (1, 'Y'): 9}
BRK = {
 (2, 'X'): [-9601109148, -8942317572, -7608305633, -6503394794, -5185864065,
            -4015634158, -3108806646, -2067921826, -1589849496, 1531948062,
            2072448179, 3208186484, 4381119427, 5851354199, 7269030693,
            8390277402, 9332614564, 9709786219],
 (2, 'Y'): [176402802, 330649406, 764952882, 1227250753, 2149465998,
            3048189112, 4322699096, 5564757994, 6801929373, 8031988371,
            8877037851, 9577761832, 9838463999, 9972727815, 9992037196],
 (1, 'X'): [-9917299785, -2259572153, 2543808026, 4437270259, 6348298970],
 (1, 'Y'): [532669786, 2504239325, 5701738806, 7966939383, 9454490138]}
# Table 2: (10^10 x, upper bound); B rows in bracket order, then P rows
TAB2 = {
 (2, 'X'): ['-0.984034048775', '-0.984375363807', '-0.984034052640', '-0.984105522615',
            '-0.984034053414', '-0.984092061375', '-0.984034037901', '-0.984364031075',
            '-0.984034038450', '-0.984033385315', '-0.984776982913', '-0.984034026898',
            '-0.984264010107', '-0.984034029391', '-0.984245044505', '-0.984034016294',
            '-0.984567630027', '-0.984033926884'],
 (2, 'Y'): ['-1.608946411646', '-1.609949505120', '-1.608960295828', '-1.609094597854',
            '-1.608960426500', '-1.608992193672', '-1.608960428486', '-1.608979509517',
            '-1.608960430478', '-1.608990548812', '-1.608960429811', '-1.609055802895',
            '-1.608960412835', '-1.609694899071', '-1.608958123669'],
 (1, 'X'): ['-2.778491531574', '-1.324666731948', '-1.324655807046', '-1.352236629957',
            '-1.324666329425'],
 (1, 'Y'): ['-1.428286151250', '-1.515602647362', '-1.428335732372', '-1.465686164672',
            '-1.428335358167']}
TAB2P = {(2, 'X'): [(-10000000000, '-0.986727371546')],
         (2, 'Y'): [(2500000000, '-1.608973617030'), (5000000000, '-1.608971701898'),
                    (7500000000, '-1.608976908052')],
         (1, 'X'): [(-10000000000, '-2.775818077526'), (-5000000000, '-1.544057819905'),
                    (5000000000, '-1.347557680876')],
         (1, 'Y'): []}
NORMTAB = {2: '.778415976284', 1: '.931985203901'}
CUT = {2: ('.77844', '-.98399', '-1.60890'), 1: ('.9321', '-1.3244', '-1.4280')}

# ------------------------------------------------------------ Q(i) arithmetic
class GQ:
    __slots__ = ('re', 'im')
    def __init__(s, re, im=Fr(0)): s.re = Fr(re); s.im = Fr(im)
    def __add__(s, o): o = o if isinstance(o, GQ) else GQ(o); return GQ(s.re+o.re, s.im+o.im)
    __radd__ = __add__
    def __neg__(s): return GQ(-s.re, -s.im)
    def __sub__(s, o): return s + (-(o if isinstance(o, GQ) else GQ(o)))
    def __mul__(s, o):
        o = o if isinstance(o, GQ) else GQ(o)
        return GQ(s.re*o.re - s.im*o.im, s.re*o.im + s.im*o.re)
    __rmul__ = __mul__
    def inv(s):
        d = s.re*s.re + s.im*s.im; return GQ(s.re/d, -s.im/d)
    def conj(s): return GQ(s.re, -s.im)
    def iszero(s): return s.re == 0 and s.im == 0
ZERO, ONE = GQ(0), GQ(1)

def trim(a):
    while len(a) > 1 and a[-1].iszero(): a = a[:-1]
    return a
def pmul(a, b):
    out = [GQ(0) for _ in range(len(a)+len(b)-1)]
    for i, x in enumerate(a):
        if x.iszero(): continue
        for j, y in enumerate(b):
            out[i+j] = out[i+j] + x*y
    return trim(out)
def padd(a, b):
    n = max(len(a), len(b))
    return trim([(a[i] if i < len(a) else ZERO) + (b[i] if i < len(b) else ZERO) for i in range(n)])
def pscale(a, c): return trim([c*x for x in a])
def pdiv_exact(a, d, tag):
    a = list(a); d = trim(d); q = [GQ(0) for _ in range(max(1, len(a)-len(d)+1))]
    li = d[-1].inv()
    for k in range(len(a)-len(d), -1, -1):
        coef = a[k+len(d)-1]*li; q[k] = coef
        for j, y in enumerate(d):
            a[k+j] = a[k+j] - coef*y
    rem_ok = all(x.iszero() for x in a)
    gate("G1", rem_ok, f"inexact division {tag}")
    return trim(q)

X1 = [ZERO, ONE]                      # x
OMX = [ONE, GQ(-1)]                   # 1 - x
OPX2 = [ONE, ZERO, ONE]               # 1 + x^2
def quad_fac(z): return [ONE + z*z, GQ(-2)*z]   # 1 - 2xz + z^2
def lin_fac(z): return [ONE, -z]                # 1 - xz

def terms(kap, which):
    """list of (z, r_z) with conjugates expanded, all exact."""
    out = []
    for re, im, a, b in EXP[(kap, which)]:
        z = GQ(Fr(re), Fr(im))
        if b is None:
            out.append((z, GQ(a*E8)))
        else:
            r = GQ(Fr(a, 2)*E8, Fr(-b, 2)*E8)
            out.append((z, r)); out.append((z.conj(), r.conj()))
    return out
def fin(kap, which): return [Fr(x)*E8 for x in FIN[(kap, which)]]

def chebU_polys(n):
    U = [[ONE], [ZERO, GQ(2)]]
    for k in range(2, n+1):
        U.append(padd(pmul([ZERO, GQ(2)], U[-1]), pscale(U[-2], GQ(-1))))
    return U
UP = chebU_polys(12)

# ------------------------------------------------------------ numerators (G1)
def build(kap):
    lp, lv = fin(kap, 'p'), fin(kap, 'v')
    tp, tv = terms(kap, 'p'), terms(kap, 'v')
    lam = LAM[kap]
    # X numerator
    QX = pmul(pmul(X1, OMX), [ONE])
    for _ in range(3-kap): QX = pmul(QX, OPX2)
    for z, _ in tp: QX = pmul(QX, quad_fac(z))
    for z, _ in tv: QX = pmul(QX, lin_fac(z))
    A = pscale(pdiv_exact(QX, X1, "QX/x"), GQ(ALPHA+2*GAMMA))
    A = padd(A, pscale(pdiv_exact(QX, OMX, "QX/(1-x)"), GQ(-2*ETA)))
    A = padd(A, pscale(pmul(X1, pdiv_exact(QX, OPX2, "QX/(1+x^2)")), GQ(-(kap+2*ALPHA+2*ETA+2*GAMMA))))
    if lam != 0:
        A = padd(A, pscale(pmul(X1, pdiv_exact(pdiv_exact(QX, OPX2, "a"), OPX2, "QX/(1+x^2)^2")), GQ(-4*lam)))
    polyt = [ZERO]
    for k, l in enumerate(lp, start=1): polyt = padd(polyt, pscale(UP[k-1], GQ(l)))
    A = padd(A, pscale(pmul(polyt, QX), GQ(-2*kap)))
    for z, r in tp:
        A = padd(A, pscale(pdiv_exact(QX, quad_fac(z), "QX/quad"), GQ(-2*kap)*r*z))
    polyh = [GQ(0) for _ in range(len(lv))]
    for k, l in enumerate(lv, start=1): polyh[k-1] = GQ(l)
    A = padd(A, pscale(pmul(polyh, QX), GQ(-1)))
    for z, r in tv:
        A = padd(A, pscale(pdiv_exact(QX, lin_fac(z), "QX/lin"), -(r*z)))
    # Y numerator
    QY = pmul(X1, OMX)
    for z, _ in tv: QY = pmul(QY, quad_fac(z))
    B = pscale(pdiv_exact(QY, X1, "QY/x"), GQ(BETA))
    B = padd(B, pscale(pdiv_exact(QY, OMX, "QY/(1-x)"), GQ(-GAMMA)))
    polytv = [ZERO]
    for k, l in enumerate(lv, start=1): polytv = padd(polytv, pscale(UP[k-1], GQ(l)))
    B = padd(B, pscale(pmul(polytv, QY), GQ(2)))
    for z, r in tv:
        B = padd(B, pscale(pdiv_exact(QY, quad_fac(z), "QY/quad"), GQ(2)*r*z))
    out = {}
    for name, P, Q in (('X', A, QX), ('Y', B, QY)):
        imz = all(c.im == 0 for c in P) and all(c.im == 0 for c in Q)
        gate("G1", imz, f"imaginary part k={kap} {name}")
        Pr = [c.re for c in P]; Qr = [c.re for c in Q]
        gate("G1", len(Pr)-1 == DEG[(kap, name)], f"degree k={kap} {name}: {len(Pr)-1}")
        out[name] = (Pr, Qr)
    return out

NUM = {}
for kap in (2, 1):
    NUM[kap] = build(kap)
    log(f"G1 k={kap}: numerator degrees X {len(NUM[kap]['X'][0])-1}, Y {len(NUM[kap]['Y'][0])-1}")

def peval(P, x):
    acc = Fr(0)
    for c in reversed(P): acc = acc*x + c
    return acc

# ------------------------------------------------------------ mp functions
mpm.mp.dps = 40
def mpc_of(z): return mpm.mpc(mpm.mpf(z.re.numerator)/z.re.denominator, mpm.mpf(z.im.numerator)/z.im.denominator)
def mpf_of(q): return mpm.mpf(q.numerator)/q.denominator
MP = {}
for kap in (2, 1):
    for w in ('p', 'v'):
        MP[(kap, w)] = ([mpf_of(l) for l in fin(kap, w)],
                        [(mpc_of(z), mpc_of(r)) for z, r in terms(kap, w)])
def chebT_vals(x, d):
    T = [mpm.mpf(1), x]
    for k in range(2, d+1): T.append(2*x*T[-1]-T[-2])
    return T
def Tfun(u, x):
    l, ex = u; T = chebT_vals(x, len(l))
    s = mpm.fsum(l[k-1]*T[k]/k for k in range(1, len(l)+1))
    s += -mpm.mpf(1)/2*mpm.fsum(r*mpm.log(1-2*x*z+z*z) for z, r in ex) if ex else 0
    return mpm.re(s) if isinstance(s, mpm.mpc) else s
def Sfun(u, x):
    l, ex = u
    s = mpm.fsum(l[k-1]*x**k/k for k in range(1, len(l)+1))
    s += -mpm.fsum(r*mpm.log(1-x*z) for z, r in ex) if ex else 0
    return mpm.re(s) if isinstance(s, mpm.mpc) else s
def norm2(u):
    l, ex = u
    s = mpm.fsum((l[k-1]**2 + 2*l[k-1]*mpm.fsum(r*z**k for z, r in ex))/k for k in range(1, len(l)+1))
    s -= mpm.fsum(r*rr*mpm.log(1-z*zz) for z, r in ex for zz, rr in ex) if ex else 0
    return mpm.re(s) if isinstance(s, mpm.mpc) else s
a_, b_, g_, e_ = (mpf_of(q) for q in (ALPHA, BETA, GAMMA, ETA))
def Xf(kap, x):
    W = (a_+2*g_)*mpm.log(abs(x)) + 2*e_*mpm.log(1-x) - (mpm.mpf(kap)/2+a_+e_+g_)*mpm.log(1+x*x)
    D = 2*g_ - 2*x*x/(1+x*x)
    return W + mpf_of(LAM[kap])*D - 2*kap*Tfun(MP[(kap, 'p')], x) - Sfun(MP[(kap, 'v')], x)
def Yf(kap, x):
    return b_*mpm.log(x) + g_*mpm.log(1-x) + 2*Tfun(MP[(kap, 'v')], x)
FUN = {'X': Xf, 'Y': Yf}

# ------------------------------------------------------------ G2 derivative check
for kap in (2, 1):
    for name, pts in (('X', ['-0.73', '-0.21', '0.17', '0.52', '0.91']),
                      ('Y', ['0.07', '0.33', '0.5', '0.68', '0.97'])):
        P, Q = NUM[kap][name]
        worst = 0
        for s in pts:
            x = Fr(s); xm = mpf_of(x)
            exact = mpf_of(peval(P, x)/peval(Q, x))
            num = mpm.diff(lambda t: FUN[name](kap, t), xm)
            rel = abs(exact-num)/max(abs(num), mpm.mpf(1))
            worst = max(worst, rel)
            gate("G2", rel < mpm.mpf('1e-25'), f"k={kap} {name} x={s} rel={mpm.nstr(rel,3)}")
        log(f"G2 k={kap} {name}: worst relative error A/Q vs d/dx = {mpm.nstr(worst,3)}")

# ------------------------------------------------------------ G3 sign variations
def variations(P, b, c, d0):
    # coefficients of (1+t)^d0 P((b+ct)/(1+t)) = sum_j A_j (b+ct)^j (1+t)^(d0-j)
    from math import comb
    acc = [Fr(0)]*(d0+1)
    for j, Aj in enumerate(P):
        if Aj == 0: continue
        bc = [Fr(comb(j, k))*b**(j-k)*c**k for k in range(j+1)]          # (b+ct)^j
        op = [comb(d0-j, k) for k in range(d0-j+1)]                        # (1+t)^(d0-j)
        for i, x in enumerate(bc):
            if x == 0: continue
            for k, y in enumerate(op):
                acc[i+k] += Aj*x*y
    sgn = [1 if a > 0 else -1 for a in acc if a != 0]
    return sum(1 for i in range(1, len(sgn)) if sgn[i] != sgn[i-1])
res = {}
for key in DIV:
    kap, name = key; P = NUM[kap][name][0]; pts = [Fr(s) for s in DIV[key]]
    v = [variations(P, pts[i], pts[i+1], DEG[key]) for i in range(len(pts)-1)]
    gate("G3", v == VAR[key], f"{key} got {v} expected {VAR[key]}")
    res[f"k{kap}{name}_variations"] = v
    log(f"G3 k={kap} {name}: sign variations {v} (paper {VAR[key]})")

# ------------------------------------------------------------ G4 brackets
TEN = Fr(10**10)
for key in BRK:
    kap, name = key; P = NUM[kap][name][0]; pts = [Fr(s) for s in DIV[key]]
    counts = [0]*(len(pts)-1); bad = 0
    for mm in BRK[key]:
        lo, hi = Fr(mm)/TEN, Fr(mm+2)/TEN
        ok = peval(P, lo)*peval(P, hi) < 0
        gate("G4", ok, f"{key} bracket {mm}")
        bad += (not ok)
        idx = [i for i in range(len(pts)-1) if pts[i] < lo and hi < pts[i+1]]
        gate("G4", len(idx) == 1, f"{key} bracket {mm} not inside one interval")
        if idx: counts[idx[0]] += 1
    gate("G4", counts == VAR[key], f"{key} bracket counts {counts}")
    log(f"G4 k={kap} {name}: {len(BRK[key])} brackets, sign-change failures {bad}, per-interval {counts}")

# ------------------------------------------------------------ G5 Table 2 values
gaps = []
vals = {}
for key in BRK:
    kap, name = key
    rows = list(zip(BRK[key], TAB2[key])) + TAB2P[key]
    vk = []
    for mm, ub in rows:
        x = mpm.mpf(mm)/10**10
        val = FUN[name](kap, x); ubm = mpm.mpf(ub)
        gate("G5", val <= ubm, f"{key} x={mm} value {mpm.nstr(val,16)} > bound {ub}")
        gap = ubm - val; gaps.append(float(gap)); vk.append(val)
        if not (gap < mpm.mpf('1.1e-12')):
            SOFT.append(f"G5 soft {key} x={mm}: gap {mpm.nstr(gap,3)}")
    vals[key] = vk
    log(f"G5 k={kap} {name}: {len(rows)} points, all <= listed bound; max computed value {mpm.nstr(max(vk),13)}")
log(f"G5 gaps (listed - computed): min {min(gaps):.3e}, max {max(gaps):.3e}")

# ------------------------------------------------------------ G6 norms
norms = {}
for kap in (2, 1):
    nv = kap*norm2(MP[(kap, 'p')]) + norm2(MP[(kap, 'v')])/2
    norms[kap] = nv; ub = mpm.mpf(NORMTAB[kap])
    gate("G6", nv <= ub, f"k={kap} norm {mpm.nstr(nv,16)} > {NORMTAB[kap]}")
    if not (ub-nv < mpm.mpf('1.1e-12')): SOFT.append(f"G6 soft k={kap}: gap {mpm.nstr(ub-nv,3)}")
    log(f"G6 k={kap}: k||p||^2+||v||^2/2 = {mpm.nstr(nv,15)} (listed {NORMTAB[kap]}, gap {mpm.nstr(ub-nv,3)})")

# ------------------------------------------------------------ G7 cutoffs and final sums
L2LO = Fr('.693146')
finals = {}
for kap in (2, 1):
    cn, cx, cy = (Fr(s) for s in CUT[kap])
    tabX = [Fr(s) for s in TAB2[(kap, 'X')]] + [Fr(s) for _, s in TAB2P[(kap, 'X')]]
    tabY = [Fr(s) for s in TAB2[(kap, 'Y')]] + [Fr(s) for _, s in TAB2P[(kap, 'Y')]]
    gate("G7", max(tabX) <= cx and max(tabY) <= cy, f"k={kap} table vs cutoff")
    gate("G7", max(vals[(kap, 'X')]) <= mpf_of(cx) and max(vals[(kap, 'Y')]) <= mpf_of(cy), f"k={kap} computed vs cutoff")
    gate("G7", norms[kap] <= mpf_of(cn), f"k={kap} norm vs cutoff")
    tot = -Fr(11, 16)*L2LO + cn + cx + cy + Fr('.000048')
    expect = Fr('-2.290939875') if kap == 2 else Fr('-2.296789875')
    gate("G7", tot == expect, f"k={kap} final {tot}")
    # sharper version from our own computed numbers and true log 2
    sharp = -mpm.mpf(11)/16*mpm.log(2) + norms[kap] + max(vals[(kap, 'X')]) + max(vals[(kap, 'Y')]) + mpm.mpf('0.000048')
    finals[kap] = {"paper": str(float(tot)), "recomputed_with_true_log2": mpm.nstr(sharp, 12)}
    log(f"G7 k={kap}: final constant {float(tot)} (paper {float(expect)}); from computed values {mpm.nstr(sharp,12)}")
lower = -Fr(8609, 4608) - Fr(2809, 4608)*Fr(693149, 10**6)
log(f"     finite-place lower constant (Prop 3.4) {float(lower):.9f}; gap to upper {float(lower)-float(Fr('-2.290939875')):.3e}")

# ------------------------------------------------------------ G8 dense-grid control (float)
def fl_terms(kap, w):
    l = [float(x) for x in fin(kap, w)]
    ex = [(complex(float(z.re), float(z.im)), complex(float(r.re), float(r.im))) for z, r in terms(kap, w)]
    return l, ex
FL = {(k, w): fl_terms(k, w) for k in (2, 1) for w in ('p', 'v')}
af, bf, gf, ef = (float(q) for q in (ALPHA, BETA, GAMMA, ETA))
def T_fl(u, x):
    l, ex = u; t0, t1 = 1.0, x; s = 0.0
    for k in range(1, len(l)+1):
        s += l[k-1]*t1/k; t0, t1 = t1, 2*x*t1-t0
    return s - 0.5*sum((r*cmath.log(1-2*x*z+z*z)).real for z, r in ex)
def S_fl(u, x):
    l, ex = u
    return sum(l[k-1]*x**k/k for k in range(1, len(l)+1)) - sum((r*cmath.log(1-x*z)).real for z, r in ex)
def X_fl(kap, x):
    W = (af+2*gf)*math.log(abs(x)) + 2*ef*math.log(1-x) - (kap/2+af+ef+gf)*math.log(1+x*x)
    D = 2*gf - 2*x*x/(1+x*x)
    return W + float(LAM[kap])*D - 2*kap*T_fl(FL[(kap, 'p')], x) - S_fl(FL[(kap, 'v')], x)
def Y_fl(kap, x):
    return bf*math.log(x) + gf*math.log(1-x) + 2*T_fl(FL[(kap, 'v')], x)
NG = 40000
grid = {}
for kap in (2, 1):
    xs = [-1 + (2*k+1)/NG for k in range(NG)]
    gx = max(X_fl(kap, x) for x in xs + [-1.0] if x != 0)
    ys = [(k+0.5)/NG for k in range(NG)]
    gy = max(Y_fl(kap, y) for y in ys)
    tabX = max(float(Fr(s)) for s in TAB2[(kap, 'X')] + [s for _, s in TAB2P[(kap, 'X')]])
    tabY = max(float(Fr(s)) for s in TAB2[(kap, 'Y')] + [s for _, s in TAB2P[(kap, 'Y')]]) if TAB2[(kap, 'Y')] else None
    gate("G8", gx <= tabX + 0.000024, f"k={kap} grid max X {gx} vs table {tabX}")
    gate("G8", gy <= tabY + 0.000024, f"k={kap} grid max Y {gy} vs table {tabY}")
    grid[kap] = {"gridmaxX": gx, "tablemaxX": tabX, "gridmaxY": gy, "tablemaxY": tabY}
    log(f"G8 k={kap}: grid max X {gx:.9f} (table max {tabX:.9f}); grid max Y {gy:.9f} (table max {tabY:.9f})")

# ------------------------------------------------------------ G9 derivative bound at bracket ends
dmax = 0.0
for key in BRK:
    kap, name = key; P, Q = NUM[kap][name]
    for mm in BRK[key]:
        x = Fr(mm)/TEN
        d = abs(float(peval(P, x)/peval(Q, x))); dmax = max(dmax, d)
        gate("G9", d < 12852, f"{key} x={mm} |deriv|={d}")
log(f"G9 max |derivative| at bracket left endpoints = {dmax:.3f} (paper bound 12852)")

total = sum(FAILS.values())
out = {"script": "cat_cert7.py", "gates": FAILS, "soft_notes": SOFT,
       "variations": res, "norms": {k: mpm.nstr(v, 16) for k, v in norms.items()},
       "table2_gap_min": min(gaps), "table2_gap_max": max(gaps),
       "finals": finals, "grid": grid, "max_abs_derivative_at_brackets": dmax,
       "GATE_FAILS": total, "seconds": round(time.time()-T0, 1)}
with open(_os.path.join(HERE, "cat_cert7_out.json"), "w") as fh:
    json.dump(out, fh, indent=1)
for s in SOFT: log(s)
log(f"gates: {FAILS}")
print(f"GATE_FAILS = {total}")
LOGF.write(f"GATE_FAILS = {total}\n"); LOGF.close()
print(f"log saved to {LOGPATH}")
