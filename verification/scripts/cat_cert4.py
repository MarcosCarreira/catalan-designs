"""
cat_cert4.py -- independent rebuild of the Section 4.4 certificate of
"Catalan's constant is irrational" (OpenAI, 24 Sep 2026), pp. 19-24.

WHAT IT CHECKS (pre-registered; expected results in brackets)
  The paper defines a fixed rational 49x48 matrix B (eq. 47) from the
  moment arrays M0, Z0 and the N=1 row polynomials P_r, D_r, and claims
  det B0, det(B0+B1), det(B0-B1) != 0 via Gaussian pivots mod 101 (Table 1).

  G1  m_i recursion == closed form 2/((i+1) c_{i/2})              [0 fails]
  G2  Y_ij (recursion of p.23) == explicit M0 formula (14)        [0 fails]
  G3  k-/k+ recursions == closed forms (19)                       [0 fails]
  G4  y_r, z_r (p.23 Chebyshev recipe) == definition (3) with N=1 [0 fails]
  G5  zeta(2) contact: coefficient of zeta(2) in every raw entry
      F_j(r), 0<=r<=48, 7<=j<=58, is exactly 0 (Prop 2.1)         [0 fails]
  G6  every denominator in B has all prime factors <= 65          [0 fails]
  G7  pivots mod 101 for sigma = 0, +1, -1 == Table 1 (144 numbers),
      swaps: none for 0 and -1; (30,31),(45,46) for +1            [0 fails]
  G8  exact rational determinants of B0, B0+B1, B0-B1 are nonzero [0 fails]
  G9  positive control vs numerical integrals of definition (6):
      M(i,j), Z(i,j) at 6 index pairs, 25 digits, with G, zeta(2) [0 fails]

RUN: python3 cat_cert4.py     (no arguments; ~10 s; single process)
OUTPUT: progress lines on screen AND in <name>_log_<YYYYMMDD_HHMMSS>.txt; cat_cert4_out.json;
        final GATE_FAILS line. All files are written next to the script.
"""
import sys, json, time
from fractions import Fraction as Fr
from math import comb, gcd
try:
    sys.set_int_max_str_digits(0)
except AttributeError:
    pass

T0 = time.time()
import os as _os
HERE = _os.path.dirname(_os.path.abspath(__file__))
STAMP = time.strftime('%Y%m%d_%H%M%S')
LOGPATH = _os.path.join(HERE, 'cat_cert4_log_' + STAMP + '.txt')
LOGF = open(LOGPATH, 'w')
def log(s):
    line = f"[{time.time()-T0:6.1f}s] {s}"
    print(line, flush=True)
    LOGF.write(line + '\n'); LOGF.flush()

FAILS = {}
def gate(name, ok, info=""):
    FAILS[name] = FAILS.get(name, 0) + (0 if ok else 1)
    if not ok:
        log(f"FAIL {name} {info}")

IMAX, JMAX = 64, 58           # array ranges used on p.23
n0, b0, q0, g0, h0, C0, L0, H0 = 48, 7, 4, 4, 2, 63, 59, 65

def c(l):                     # c_l = 4^-l binom(2l,l); zero unless l integer >= 0
    return Fr(comb(2*l, l), 4**l)

# ---------- moments m_i  (G1)
m = [Fr(0)]*(IMAX+1)
m[0] = Fr(2)
for i in range(2, IMAX+1, 2):
    m[i] = i*m[i-2]/(i+1)
for i in range(IMAX+1):
    ref = Fr(2)/((i+1)*c(i//2)) if i % 2 == 0 else Fr(0)
    gate("G1", m[i] == ref, f"i={i}")
def mm(i):                    # m_{-1} = 0 convention
    return m[i] if 0 <= i <= IMAX else (Fr(0) if i == -1 else None)
log("G1 moments done")

# ---------- boundary arrays k-, k+  (G3)
km = [Fr(0)]*(IMAX+1); km[1] = Fr(2)
for d in range(2, IMAX+1):
    km[d] = ((d-1)*km[d-2] + m[d-2] + m[d-1])/d
kp = [Fr(0)]*(JMAX+1)
for d in range(2, JMAX+1):
    kp[d] = ((d-2)*kp[d-2] + Fr(2, d-1))/(d-1)
def Hs(z):                    # H*_z of (18)
    return c(z//2) if z % 2 == 0 else Fr(1, z)/c((z-1)//2)
for u in range(IMAX+1):
    if u % 2 == 0:
        ref = Hs(u)*sum((Fr(2, z*z)/(Hs(z)**2) for z in range(2, u+1, 2)), Fr(0))
    else:
        ref = Hs(u)*sum((Fr(2, z) for z in range(1, u+1, 2)), Fr(0))
    gate("G3", km[u] == ref, f"k-_{u}")
for u in range(0, JMAX):
    ref = Hs(u)*sum((Fr(2, z*z)/Hs(z) for z in range(1, u+1) if (z-u) % 2 == 0), Fr(0))
    gate("G3", kp[u+1] == ref, f"k+_{u+1}")
log("G3 boundary arrays done")

# ---------- Y = M0 array (p.23) and explicit (14)  (G2)
Y = [[Fr(0)]*(JMAX+1) for _ in range(IMAX+1)]
for i in range(IMAX+1): Y[i][0] = km[i]
for j in range(JMAX+1): Y[0][j] = kp[j]
for i in range(1, IMAX+1):
    for j in range(1, JMAX+1):
        Y[i][j] = Y[i-1][j-1] - m[i-1]/j
for i in range(IMAX+1):
    for j in range(JMAX+1):
        if i >= j:
            ref = km[i-j] - sum((mm(i-j+k-1)/k for k in range(1, j+1)), Fr(0))
        else:
            ref = kp[j-i] - sum((mm(k-(j-i)-1)/k for k in range(j-i+1, j+1)), Fr(0))
        gate("G2", Y[i][j] == ref, f"({i},{j})")
log("G2 M0 array done")

# ---------- Z0 and J = (3/2) Z0
B1 = [Fr(0)]*(IMAX+1); B2 = [Fr(0)]*(IMAX+1)
for i in range(1, IMAX+1):
    B1[i] = B1[i-1] + Fr(1, i); B2[i] = B2[i-1] + Fr(1, i*i)
def Z0(i, j):
    return -B2[i] if i == j else (B1[i]-B1[j])/(i-j)
J = [[Fr(3, 2)*Z0(i, j) for j in range(JMAX+1)] for i in range(IMAX+1)]

# ---------- row polynomials (coefficient lists, degree <= 64)
DEG = 64
def pmul(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i+j] += x*y
    return out
def shift(a, k):              # multiply by t^k (k may be negative if exact)
    if k >= 0: return [0]*k + a
    assert all(x == 0 for x in a[:-k]); return a[-k:]
def pad(a, n=DEG+1):
    a = a + [0]*(n-len(a)); assert all(x == 0 for x in a[n:]); return a[:n]
def mono(k): return [0]*k + [1]
E = [mono(62), mono(61)]; I = [[0], mono(62)]
for mi in range(2, 45):
    for S in (E, I):
        a = shift(S[mi-1], -1); a = [2*x for x in a]
        b = S[mi-2]; L = max(len(a), len(b))
        S.append([ (a[k] if k < len(a) else 0) - (b[k] if k < len(b) else 0) for k in range(L)])
oneminus2 = [1, -2, 1]
yr, zr = [], []
for r in range(49):
    u = abs(r-4); sg = (r > 4) - (r < 4)
    yr.append(pad(pmul(oneminus2, E[u])))
    zr.append(pad([sg*x for x in pmul(oneminus2, I[u])]))
# (G4) independent check against (3): Chebyshev T_d, U_{d-1} via w-Laurent
def cheb(dmax):
    T = [[1], [0, 1]]; U = [[1], [0, 2]]
    for k in range(2, dmax+1):
        T.append([2*x for x in [0]+T[-1]]); T[-1] = [ (T[-1][i] if i < len(T[-1]) else 0) - (T[-3][i] if i < len(T[-3]) else 0) for i in range(len(T[-1]))]
        U.append([2*x for x in [0]+U[-1]]); U[-1] = [ (U[-1][i] if i < len(U[-1]) else 0) - (U[-3][i] if i < len(U[-3]) else 0) for i in range(len(U[-1]))]
    return T, U
Tc, Uc = cheb(50)
for r in range(49):
    d = abs(r-g0); sg = (r > g0) - (r < g0)
    # (1-t)^h t^{C-1} T_d(1/t): coefficient of x^k in T_d goes to t^{C-1-k}
    P = [0]*(C0+1)
    for k, x in enumerate(Tc[d]): P[C0-1-k] += x
    D = [0]*(C0+1)
    if d >= 1:
        for k, x in enumerate(Uc[d-1]): D[C0-1-k] += sg*x
    gate("G4", pad(pmul(oneminus2, P)) == yr[r], f"P r={r}")
    gate("G4", pad(pmul(oneminus2, D)) == zr[r], f"D r={r}")
log("G4 row polynomials done")

# ---------- zeta(2) contact (G5): sum_i [t^i]P c_{(j-i-1)/2} - [t^j]D == 0
for r in range(49):
    for j in range(b0, L0):
        s = Fr(0)
        for i in range(0, j):
            if (j-i-1) % 2 == 0 and yr[r][i]:
                s += yr[r][i]*c((j-i-1)//2)
        s -= zr[r][j]
        gate("G5", s == 0, f"r={r} j={j} value={s}")
log("G5 zeta(2) cancellation done")

# ---------- v_r(j) and B (eq. 47)
B = []
for r in range(49):
    v = [sum((yr[r][i]*Y[i][j] - zr[r][i]*J[i][j] for i in range(DEG+1)), Fr(0)) for j in range(7, 59)]
    for _ in range(4):
        v = [v[k]-v[k+1] for k in range(len(v)-1)]
    assert len(v) == 48
    B.append(v)
log("B (49x48) built")

def primefactors(nn):
    nn = abs(nn); ps = set(); p = 2
    while p*p <= nn:
        while nn % p == 0: ps.add(p); nn //= p
        p += 1
    if nn > 1: ps.add(nn)
    return ps
maxp = 0
for row in B:
    for x in row:
        ps = primefactors(x.denominator); mp_ = max(ps) if ps else 1
        maxp = max(maxp, mp_)
        gate("G6", mp_ <= 65, f"den prime {mp_}")
log(f"G6 largest denominator prime in B = {maxp}")

# ---------- pivots mod 101 (G7)
P101 = 101
TABLE1 = {
 0: [60,68,62,79,47,32,69,57,30,72,35,66,43,20,85,48,
     88,4,77,54,60,79,26,68,83,39,40,65,1,68,78,24,
     15,98,32,22,94,9,99,10,15,75,4,2,25,53,90,79],
 1: [38,51,90,70,5,21,88,55,45,20,35,41,77,10,18,25,
     76,14,38,72,6,66,56,35,83,58,56,11,17,20,30,24,
     28,11,25,70,79,99,66,38,4,41,63,91,63,17,90,98],
 -1:[82,92,26,87,21,74,87,88,88,3,14,23,38,58,36,20,
     26,33,94,74,78,45,93,86,73,66,45,30,61,3,88,27,
     20,58,69,48,78,39,48,1,66,18,86,93,52,92,39,77]}
SWAPS = {0: [], 1: [(30, 31), (45, 46)], -1: []}
def modp(x, p): return x.numerator*pow(x.denominator, -1, p) % p
res7 = {}
for sg in (0, 1, -1):
    R = [[modp(B[i][k] + sg*B[i+1][k], P101) for k in range(48)] for i in range(48)]
    piv, swaps = [], []
    for i in range(48):
        if R[i][i] == 0:
            jj = next((j for j in range(i+1, 48) if R[j][i] != 0), None)
            if jj is None: piv.append(0); break
            R[i], R[jj] = R[jj], R[i]; swaps.append((i, jj))
        di = R[i][i]; piv.append(di); inv = pow(di, -1, P101)
        for j in range(i+1, 48):
            f = R[j][i]*inv % P101
            if f: R[j] = [(R[j][k]-f*R[i][k]) % P101 for k in range(48)]
    ok = (piv == TABLE1[sg]) and (swaps == SWAPS[sg])
    gate("G7", ok, f"sigma={sg} pivots={piv} swaps={swaps}")
    nmatch = sum(1 for a, b in zip(piv, TABLE1[sg]) if a == b)
    res7[str(sg)] = {"pivots_match": nmatch, "swaps": swaps}
    log(f"G7 sigma={sg:+d}: {nmatch}/48 pivots match Table 1, swaps {swaps}")

# ---------- exact determinants over Q (G8): Bareiss on integer-scaled rows
def det_exact(Mq):
    n = len(Mq); A = []; scale = Fr(1)
    for row in Mq:
        l = 1
        for x in row: l = l*x.denominator//gcd(l, x.denominator)
        A.append([int(x*l) for x in row]); scale /= l
    sign = 1; prev = 1
    for k in range(n-1):
        if A[k][k] == 0:
            sw = next((i for i in range(k+1, n) if A[i][k] != 0), None)
            if sw is None: return Fr(0)
            A[k], A[sw] = A[sw], A[k]; sign = -sign
        for i in range(k+1, n):
            for j in range(k+1, n):
                A[i][j] = (A[i][j]*A[k][k] - A[i][k]*A[k][j])//prev
        prev = A[k][k]
    return sign*A[n-1][n-1]*scale
res8 = {}
for sg in (0, 1, -1):
    d = det_exact([[B[i][k] + sg*B[i+1][k] for k in range(48)] for i in range(48)])
    gate("G8", d != 0, f"sigma={sg}")
    res8[str(sg)] = {"num_digits": len(str(abs(d.numerator))), "den_digits": len(str(d.denominator)),
                     "mod101": modp(d, P101) if d.denominator % P101 else None}
    log(f"G8 sigma={sg:+d}: det != 0, numerator {res8[str(sg)]['num_digits']} digits, "
        f"denominator {res8[str(sg)]['den_digits']} digits")

# ---------- numerical positive control (G9)
import mpmath as mpm
mpm.mp.dps = 30
Gc, z2 = mpm.catalan, mpm.zeta(2)
def gj(t, j, l1mt=None):      # int_0^1 s^j/(1-ts) ds ; l1mt = log(1-t) if supplied
    if abs(t) < mpm.mpf("0.5"):
        return mpm.fsum(t**k/(j+k+1) for k in range(120))
    L = mpm.log(1-t) if l1mt is None else l1mt
    return -(L + sum(t**k/k for k in range(1, j+1)))/t**(j+1)
def Mnum(i, j):
    # t = sin(theta): negative half directly; positive half via t = cos(d),
    # 1 - t = 2 sin(d/2)^2, so log(1-t) never rounds to log(0) at t -> 1
    neg = mpm.quad(lambda th: mpm.sin(th)**i*abs(mpm.sin(th))*gj(mpm.sin(th), j), [-mpm.pi/2, 0])
    def fpos(d):
        t = mpm.cos(d)
        return t**i*t*gj(t, j, mpm.log(2*mpm.sin(d/2)**2))
    return neg + mpm.quad(fpos, [0, mpm.pi/2])
def Znum(i, j):
    return mpm.quad(lambda t: t**i*gj(t, j), [0, 1])
def cc(z):                    # c_z, zero unless z integer >= 0
    return mpm.mpf(comb(2*z, z))/4**z if z == int(z) and z >= 0 else mpm.mpf(0)
pos = []
for (i, j) in [(0, 0), (3, 0), (0, 3), (2, 5), (5, 2), (4, 4)]:
    ex_M = mpm.mpf(Y[i][j].numerator)/Y[i][j].denominator \
        + 4*Gc*(cc((i-j)//2) if (i-j) % 2 == 0 and i >= j else 0) \
        + mpm.mpf(3)/2*z2*(cc((j-i-1)//2) if (j-i-1) % 2 == 0 and j > i else 0)
    zz = Z0(i, j); ex_Z = mpm.mpf(zz.numerator)/zz.denominator + (z2 if i == j else 0)
    eM = abs(Mnum(i, j)-ex_M); eZ = abs(Znum(i, j)-ex_Z)
    gate("G9", eM < 1e-20 and eZ < 1e-20, f"({i},{j}) errM={eM} errZ={eZ}")
    pos.append({"ij": [i, j], "errM": mpm.nstr(eM, 3), "errZ": mpm.nstr(eZ, 3)})
    log(f"G9 ({i},{j}): |M_num-M_formula| = {mpm.nstr(eM,3)}, |Z_num-Z_formula| = {mpm.nstr(eZ,3)}")

total = sum(FAILS.values())
out = {"script": "cat_cert4.py", "gates": FAILS, "G6_max_den_prime": maxp,
       "G7": res7, "G8": res8, "G9": pos, "GATE_FAILS": total,
       "seconds": round(time.time()-T0, 1)}
with open(_os.path.join(HERE, "cat_cert4_out.json"), "w") as fh:
    json.dump(out, fh, indent=1)
log(f"gates: {FAILS}")
print(f"GATE_FAILS = {total}")
LOGF.write(f"GATE_FAILS = {total}\n"); LOGF.close()
print(f"log saved to {LOGPATH}")
