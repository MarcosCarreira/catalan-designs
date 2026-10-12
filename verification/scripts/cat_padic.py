"""
cat_padic.py -- computational checks of the finite-place (lower-bound) side of
"Catalan's constant is irrational" (OpenAI, 24 Sep 2026), Secs. 3.1-3.3, plus an
end-to-end look at the actual determinants Delta_N for N = 1, 2.

WHAT IT CHECKS (pre-registered; expected results in brackets)
  A  Lemma 3.1 by brute force: for p = 7, 11, 13, 17, 19, 23 and the largest
     allowed range H = (p^2-1)//4, every pair 0 <= i, j < H:
       p^2 M0(i,j) is p-integral and == E_d M0(i',j') (mod p);
       p^2 Z0(i,j) is p-integral and == Z0(i//p, j//p) or 0 (mod p);
     also (34) and both congruences (35) for every u < H.          [0 fails]
  B  Lemma 3.3 counting identity R0 = K_+ + (J-b)_+ - B equals the stated R,
     for N = 1..60 and every prime H/2 < p <= H.                    [0 fails]
  C  END-TO-END per-prime bound. Delta_N(G) is computed exactly (as a function
     of a hypothetical rational value of G) for N = 1 (n = 48) and N = 2 (n = 96),
     with G = 0 and G = 1/3. For every prime 2 sqrt(H) < p <= H, the p-adic
     valuation is computed by elimination over Z/p^K (full pivoting by
     valuation, exact) and compared with the paper's per-prime loss bound:
       p <= H/2:  -v_p(Delta) <= 2n - (d0 - q)_+          (eq. 38)
       p >  H/2:  -v_p(Delta) <= min(2n, n+R, 2R+S)        (eq. 40)  [0 fails]
     positive control: primes H < p < H + 60 give v_p >= 0           [0 fails]
     report: slack per prime; v_p for small primes p <= 2 sqrt(H) (no claim).
  D  REAL PLACE at finite N (report only, no gate): Delta_N with the true G,
     computed at two working precisions (agreement to 30 digits is the control
     [0 fails]); prints log|Delta_N|/n^2 - (1/2) log 2. The paper's limsup bound is
     -2.290939875 (Prop. 7.1) and the liminf under rationality is > -2.29084
     (Prop. 3.4); at N = 1, 2 finite-size terms O(log n / n) dominate, so these
     numbers show the trend only.

RUN: python3 cat_padic.py     (no arguments; single process; needs mpmath)
TIME (record run on Marcos's Mac, 2026-10-06): A 0.3 s, B < 0.1 s, C N=1 15 s, C N=2 11 min, D 5 s. Total ~11.5 min.
OUTPUT: progress lines on screen AND in <name>_log_<YYYYMMDD_HHMMSS>.txt; cat_padic_out.json;
        final GATE_FAILS line. All files are written next to the script.
"""
import sys, json, time, math
from fractions import Fraction as Fr
from math import comb
try:
    sys.set_int_max_str_digits(0)
except AttributeError:
    pass
import mpmath as mpm

NS = [1, 2]                      # determinant sizes n = 48 N
GVALS = [Fr(0), Fr(1, 3)]        # hypothetical rational values of G
PA = [7, 11, 13, 17, 19, 23]     # primes for Part A
DPS = {1: (1000, 1500), 2: (1500, 2500)}   # Part D working precisions (agreement is the control)

T0 = time.time()
import os as _os
HERE = _os.path.dirname(_os.path.abspath(__file__))
STAMP = time.strftime('%Y%m%d_%H%M%S')
LOGPATH = _os.path.join(HERE, 'cat_padic_log_' + STAMP + '.txt')
LOGF = open(LOGPATH, 'w')
def log(s):
    line = f"[{time.time()-T0:7.1f}s] {s}"
    print(line, flush=True)
    LOGF.write(line + '\n'); LOGF.flush()
FAILS = {}
def gate(name, ok, info=""):
    FAILS[name] = FAILS.get(name, 0) + (0 if ok else 1)
    if not ok: log(f"FAIL {name} {info}")

def c(l): return Fr(comb(2*l, l), 4**l)
def vp(x, p):
    if x == 0: return 10**9
    v = 0; a, b = x.numerator, x.denominator
    while a % p == 0: a //= p; v += 1
    while b % p == 0: b //= p; v -= 1
    return v
def primes_upto(m):
    s = [True]*(m+1); s[0] = s[1] = False
    for i in range(2, int(m**0.5)+1):
        if s[i]:
            for j in range(i*i, m+1, i): s[j] = False
    return [i for i in range(m+1) if s[i]]

def arrays(Imax, Jmax):
    """m_i (i <= Imax), k-_d (d <= Imax), k+_d (d <= Jmax+1), M0[i][j], Z0 via B1, B2."""
    top = max(Imax, Jmax) + 2
    m = [Fr(0)]*(top+1); m[0] = Fr(2)
    for i in range(2, top+1, 2): m[i] = i*m[i-2]/(i+1)
    km = [Fr(0)]*(top+1); km[1] = Fr(2)
    for d in range(2, top+1): km[d] = ((d-1)*km[d-2] + m[d-2] + m[d-1])/d
    kp = [Fr(0)]*(top+1)
    for d in range(2, top+1): kp[d] = ((d-2)*kp[d-2] + Fr(2, d-1))/(d-1)
    M = [[Fr(0)]*(Jmax+1) for _ in range(Imax+1)]
    for i in range(Imax+1): M[i][0] = km[i]
    for j in range(Jmax+1): M[0][j] = kp[j]
    for i in range(1, Imax+1):
        Mi, Mp, mi = M[i], M[i-1], m[i-1]
        for j in range(1, Jmax+1): Mi[j] = Mp[j-1] - mi/j
    B1 = [Fr(0)]*(top+1); B2 = [Fr(0)]*(top+1)
    for i in range(1, top+1): B1[i] = B1[i-1] + Fr(1, i); B2[i] = B2[i-1] + Fr(1, i*i)
    def Z0(i, j): return -B2[i] if i == j else (B1[i]-B1[j])/(i-j)
    return m, km, kp, M, Z0

# =============================================================== Part A
log("Part A: Lemma 3.1 brute force")
resA = {}
for p in PA:
    H = (p*p - 1)//4
    m, km, kp, M, Z0 = arrays(H, H)
    Hr = H//p + 2
    def M0red(i_, j_):
        return kp[j_+1] if i_ == -1 else M[i_][j_]
    def E(d):
        return c(d//2) if (0 <= d < p and d % 2 == 0) else Fr(0)
    bad = 0
    for i in range(H):
        for j in range(H):
            l = j % p; d = (j-i-1) % p
            ip = (i+1+d-l)//p - 1; jp = j//p
            assert (i+1+d-l) % p == 0 and ip >= -1
            x = p*p*M[i][j]; y = E(d)*M0red(ip, jp)
            ok = vp(x, p) >= 0 and vp(x-y, p) >= 1
            z = p*p*Z0(i, j); w = Z0(i//p, j//p) if (i-j) % p == 0 else Fr(0)
            ok = ok and vp(z, p) >= 0 and vp(z-w, p) >= 1
            if not ok:
                bad += 1
                if bad <= 5: log(f"FAIL A p={p} (i,j)=({i},{j})")
    eps = (-1)**((p-1)//2)
    for u in range(H):
        P_, r = divmod(u, p)
        # (35)
        ok = vp(p*p*km[u] - E(p-1-r)*km[P_], p) >= 1 and vp(p*p*km[u], p) >= 0
        ok = ok and vp(p*p*kp[u+1] - E(r)*kp[P_+1], p) >= 1 and vp(p*p*kp[u+1], p) >= 0
        # (34) with z = u
        z = u; Pz, rz = divmod(z, p)
        lhs = p*(m[z-1] if z >= 1 else Fr(0)); rhs = E(p-1-rz)*(m[Pz-1] if Pz >= 1 else Fr(0))
        ok = ok and vp(lhs - rhs, p) >= 1 and vp(lhs, p) >= 0
        if not ok:
            bad += 1
            if bad <= 5: log(f"FAIL A p={p} u={u} (34)/(35)")
    FAILS["A"] = FAILS.get("A", 0) + bad
    resA[p] = {"H": H, "pairs": H*H, "fails": bad}
    log(f"  p={p:2d}: H={H:3d}, {H*H} pairs (i,j) + {H} values u, failures {bad}")

# =============================================================== Part B
log("Part B: Lemma 3.3 counting identity")
pos = lambda y: max(y, 0)
nb = 0
for N in range(1, 61):
    n, a, b, q, g, h = 48*N, 11*N, 7*N, 4*N, 4*N, 2*N
    L, C, H, A = 59*N, 63*N, 65*N, 19*N
    for p in primes_upto(H):
        if 2*p <= H: continue
        K, J = L-p, H-p
        Bp = pos(min(A, K) - b)
        R0 = pos(K) + pos(J-b) - Bp
        R = pos(K) + pos(K-A) + pos(J - max(b, K))
        gate("B", R0 == R, f"N={N} p={p} R0={R0} R={R}"); nb += 1
log(f"  {nb} (N,p) cases checked")

# =============================================================== Part C
def build(N):
    n, a, b, q, g, h = 48*N, 11*N, 7*N, 4*N, 4*N, 2*N
    L, C, H, A = 59*N, 63*N, 65*N, 19*N
    m, km, kp, M, Z0 = arrays(H, L)
    # Chebyshev T_d, U_{d-1} as coefficient lists in x
    T = [[1], [0, 1]]; U = [[1], [0, 2]]
    for k in range(2, n+2):
        for S in (T, U):
            a1 = [0] + [2*v for v in S[-1]]; a0 = S[-2]
            S.append([(a1[i] if i < len(a1) else 0) - (a0[i] if i < len(a0) else 0) for i in range(len(a1))])
    binh = [(-1)**v*comb(h, v) for v in range(h+1)]
    def rowpoly(coeffs_x, sign):      # (1-t)^h t^(C-1) Poly(1/t)
        base = {}
        for k, v in enumerate(coeffs_x):
            if v: base[C-1-k] = sign*v
        out = {}
        for e, v in base.items():
            for s, bv in enumerate(binh):
                out[e+s] = out.get(e+s, 0) + v*bv
        return {e: v for e, v in out.items() if v}
    Pr, Dr = [], []
    for r in range(n):
        d = abs(r-g); sg = (r > g) - (r < g)
        Pr.append(rowpoly(T[d], 1))
        Dr.append(rowpoly(U[d-1], sg) if d >= 1 else {})
    # raw columns j in [b, L): rational part and G-coefficient
    Zc = {}
    def Z0c(i, j):
        key = (i, j)
        if key not in Zc: Zc[key] = Z0(i, j)
        return Zc[key]
    raw_rat, raw_G = [], []
    for r in range(n):
        rr, rg = [], []
        for j in range(b, L):
            s = Fr(0); sG = Fr(0)
            for i, v in Pr[r].items():
                s += v*M[i][j]
                if (i-j) % 2 == 0 and i >= j: sG += v*c((i-j)//2)
            for i, v in Dr[r].items():
                s -= Fr(3, 2)*v*Z0c(i, j)
            rr.append(s); rg.append(4*sG)
        raw_rat.append(rr); raw_G.append(rg)
    # filter: column k = sum_j (-1)^(j-b-k) C(q, j-b-k) raw_j
    Frat = [[sum(((-1)**t)*comb(q, t)*raw_rat[r][k+t] for t in range(q+1)) for k in range(n)] for r in range(n)]
    FG = [[sum(((-1)**t)*comb(q, t)*raw_G[r][k+t] for t in range(q+1)) for k in range(n)] for r in range(n)]
    return dict(n=n, a=a, b=b, q=q, g=g, h=h, L=L, C=C, H=H, A=A, Frat=Frat, FG=FG)

def vp_det(Mq, p):
    """exact p-adic valuation of det of a rational matrix (full pivoting over Z/p^K)."""
    n = len(Mq)
    s = -min(min(vp(x, p) for x in row) for row in Mq)
    s = max(s, 0)
    K = n*s + 3*n + 20; MOD = p**K
    A = []
    for row in Mq:
        rr = []
        for x in row:
            if x == 0: rr.append(0); continue
            num, den = x.numerator*(p**s), x.denominator
            while den % p == 0:
                den //= p; num //= p   # num divisible by p here since x*p^s is p-integral
            rr.append(num*pow(den, -1, MOD) % MOD)
        A.append(rr)
    def val(x):
        if x == 0: return K
        v = 0
        while x % p == 0: x //= p; v += 1
        return v
    tot = 0
    rows = list(range(n)); cols = list(range(n))
    for step in range(n):
        best = (K+1, None, None)
        for i in range(step, n):
            Ai = A[i]
            for j in range(step, n):
                x = Ai[j]
                if x:
                    v = val(x)
                    if v < best[0]:
                        best = (v, i, j)
                        if v == 0: break
            if best[0] == 0: break
        v, bi, bj = best
        if bi is None:                       # remaining block is 0 mod p^K
            return tot + (n-step)*K - n*s, True
        A[step], A[bi] = A[bi], A[step]
        for row in A: row[step], row[bj] = row[bj], row[step]
        tot += v
        piv = A[step][step]; u = piv // p**v; uinv = pow(u, -1, MOD)
        prow = A[step]
        for i in range(step+1, n):
            x = A[i][step]
            if x:
                f = (x // p**v) * uinv % MOD
                Ai = A[i]
                A[i] = [(Ai[k] - f*prow[k]) % MOD for k in range(n)]
    return tot - n*s, False

resC = {}
DATA = {}
for N in NS:
    log(f"Part C: building Delta_{N} exactly (n = {48*N})")
    D = build(N); DATA[N] = D
    n, q, b, L, H, A = D['n'], D['q'], D['b'], D['L'], D['H'], D['A']
    log(f"  built; primes 2*sqrt(H) < p <= H with H = {H}")
    for G in GVALS:
        Mq = [[D['Frat'][r][k] + G*D['FG'][r][k] for k in range(n)] for r in range(n)]
        rows = []
        for p in primes_upto(H + 60):
            if p == 2: continue
            if G.denominator % p == 0: continue
            v, lb = vp_det(Mq, p)
            if p*p <= 4*H:
                rows.append({"p": p, "vp": v, "claim": None}); continue
            if p > H:
                gate("C+", v >= 0, f"N={N} G={G} p={p} vp={v}")
                rows.append({"p": p, "vp": v, "claim": ">=0"}); continue
            if 2*p <= H:
                d0 = pos(min(p, L-p, A) - max(b, H-2*p))
                Dp = 2*n - pos(d0 - q); which = "38"
            else:
                K_, J_ = L-p, H-p
                Bp = pos(min(A, K_) - b)
                R = pos(K_) + pos(K_-A) + pos(J_ - max(b, K_))
                S = Bp + pos(min(b, J_) - max(0, K_))
                Dp = min(2*n, n+R, 2*R+S); which = "40"
            ok = -v <= Dp
            gate("C", ok, f"N={N} G={G} p={p} vp={v} bound -{Dp}")
            rows.append({"p": p, "vp": v, "bound": -Dp, "eq": which, "slack": v + Dp})
        resC[f"N{N}_G{G}"] = rows
        claimed = [x for x in rows if x.get("bound") is not None]
        log(f"  G={str(G):>3}: {len(claimed)} primes in range, all v_p >= bound: "
            f"{all(x['vp'] >= x['bound'] for x in claimed)}; slack min {min(x['slack'] for x in claimed)}, "
            f"max {max(x['slack'] for x in claimed)}; primes > H: v_p >= 0 for all: "
            f"{all(x['vp'] >= 0 for x in rows if x.get('claim') == '>=0')}")
        log("    p : v_p (bound)  " + "  ".join(f"{x['p']}:{x['vp']}({x['bound']})" for x in claimed))

# =============================================================== Part D
resD = {}
for N in NS:
    D = DATA[N]; n = D['n']
    outs = []
    for dps in DPS[N]:
        mpm.mp.dps = dps
        Gm = mpm.catalan
        Mm = mpm.matrix(n, n)
        for r in range(n):
            for k in range(n):
                x, y = D['Frat'][r][k], D['FG'][r][k]
                Mm[r, k] = mpm.mpf(x.numerator)/x.denominator + Gm*(mpm.mpf(y.numerator)/y.denominator)
        dt = mpm.det(Mm)
        outs.append(dt)
        log(f"Part D: N={N}, dps={dps}: log10|Delta| = {mpm.nstr(mpm.log10(abs(dt)), 12)}")
    mpm.mp.dps = 60
    rel = abs(outs[0]-outs[1])/abs(outs[1])
    gate("D", rel < mpm.mpf('1e-30'), f"N={N} precision disagreement {mpm.nstr(rel,3)}")
    Lval = mpm.log(abs(outs[1]))/n**2 - mpm.log(2)/2
    resD[N] = {"log10_abs_Delta": mpm.nstr(mpm.log10(abs(outs[1])), 15),
               "normalized": mpm.nstr(Lval, 12), "rel_diff_two_precisions": mpm.nstr(rel, 3),
               "sign": int(mpm.sign(outs[1]))}
    log(f"Part D: N={N}: log|Delta_N|/n^2 - log2/2 = {mpm.nstr(Lval, 10)}  "
        f"(paper: limsup <= -2.290939875; two precisions agree to {mpm.nstr(rel,3)})")

total = sum(FAILS.values())
out = {"script": "cat_padic.py", "gates": FAILS, "A": resA, "C": resC, "D": resD,
       "GATE_FAILS": total, "seconds": round(time.time()-T0, 1)}
with open(_os.path.join(HERE, "cat_padic_out.json"), "w") as fh: json.dump(out, fh, indent=1, default=str)
log(f"gates: {FAILS}")
print(f"GATE_FAILS = {total}")
LOGF.write(f"GATE_FAILS = {total}\n"); LOGF.close()
print(f"log saved to {LOGPATH}")
