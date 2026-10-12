"""full_cat.py (scout) -- measured balance of the preprint's own determinant for G, same normalisation as b4_balance.py:
lambda = log|Delta|/n^2 (true G), varpi = -log den(Delta with G -> 1/3)/n^2, gap = lambda - varpi; n = matrix size.
Construction copied from cat_extend.py build(N), generalised to integer parameters (n, b, q, g, h).
Usage: python3 full_cat.py n,b,q,g,h [n,b,q,g,h ...]   or   python3 full_cat.py ratio n1,n2,...  (preprint ratios 48:7:4:4:2)"""
import sys, math, time
from fractions import Fraction as Fr
from math import comb
import mpmath as mpm
sys.path.insert(0, '.')
try: sys.set_int_max_str_digits(0)
except AttributeError: pass
T0 = time.time()
def c(l): return Fr(comb(2*l, l), 4**l)

def arrays(Imax, Jmax):
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

def build(n, b, q, g, h):
    a = b+q; L = n+a; C = L+g; H = C+h
    m, km, kp, M, Z0 = arrays(H, L)
    T = [[1], [0, 1]]; U = [[1], [0, 2]]
    for k in range(2, n+2):
        for S in (T, U):
            a1 = [0] + [2*v for v in S[-1]]; a0 = S[-2]
            S.append([(a1[i] if i < len(a1) else 0) - (a0[i] if i < len(a0) else 0) for i in range(len(a1))])
    binh = [(-1)**v*comb(h, v) for v in range(h+1)]
    def rowpoly(coeffs_x, sign):
        base = {C-1-k: sign*v for k, v in enumerate(coeffs_x) if v}
        out = {}
        for e, v in base.items():
            for s, bv in enumerate(binh): out[e+s] = out.get(e+s, 0) + v*bv
        return {e: v for e, v in out.items() if v}
    Pr, Dr = [], []
    for r in range(n):
        d = abs(r-g); sg = (r > g) - (r < g)
        Pr.append(rowpoly(T[d], 1)); Dr.append(rowpoly(U[d-1], sg) if d >= 1 else {})
    raw_rat, raw_G = [], []
    for r in range(n):
        rr, rg = [], []
        for j in range(b, L):
            s = Fr(0); sG = Fr(0)
            for i, v in Pr[r].items():
                s += v*M[i][j]
                if (i-j) % 2 == 0 and i >= j: sG += v*c((i-j)//2)
            for i, v in Dr[r].items(): s -= Fr(3, 2)*v*Z0(i, j)
            rr.append(s); rg.append(4*sG)
        raw_rat.append(rr); raw_G.append(rg)
    Frat = [[sum(((-1)**t)*comb(q, t)*raw_rat[r][k+t] for t in range(q+1)) for k in range(n)] for r in range(n)]
    FG = [[sum(((-1)**t)*comb(q, t)*raw_G[r][k+t] for t in range(q+1)) for k in range(n)] for r in range(n)]
    return Frat, FG

def det_exact(Mx):
    n = len(Mx); scale = 1; A = []
    for row in Mx:
        Lc = 1
        for x in row: Lc = Lc*x.denominator//math.gcd(Lc, x.denominator)
        scale *= Lc; A.append([int(x*Lc) for x in row])
    sign = 1; prev = 1
    for cc in range(n-1):
        p = next((r for r in range(cc, n) if A[r][cc] != 0), None)
        if p is None: return Fr(0)
        if p != cc: A[cc], A[p] = A[p], A[cc]; sign = -sign
        for r in range(cc+1, n):
            for k in range(cc+1, n): A[r][k] = (A[r][k]*A[cc][cc] - A[r][cc]*A[cc][k])//prev
            A[r][cc] = 0
        prev = A[cc][cc]
    return Fr(sign*A[n-1][n-1], scale)

def measure(n, b, q, g, h):
    t = time.time()
    Frat, FG = build(n, b, q, g, h)
    tb = time.time()-t
    Dh = det_exact([[Frat[r][k] + Fr(1, 3)*FG[r][k] for k in range(n)] for r in range(n)])
    out = []
    for dps in (60 + 4*n*n, 110 + 4*n*n):
        mpm.mp.dps = dps
        G = +mpm.catalan
        Mm = mpm.matrix(n, n)
        for r in range(n):
            for k in range(n):
                x, y = Frat[r][k], FG[r][k]
                Mm[r, k] = mpm.mpf(x.numerator)/x.denominator + G*mpm.mpf(y.numerator)/y.denominator
        out.append(mpm.det(Mm))
    rel = abs(out[0]-out[1])/abs(out[1])
    lam = float(mpm.log(abs(out[1]))/n**2); vp = -math.log(Dh.denominator)/n**2
    return lam, vp, lam-vp, float(rel), tb, time.time()-t

if __name__ == "__main__":
    if sys.argv[1] == 'ratio':
        designs = [(n, round(7*n/48), round(4*n/48), round(4*n/48), round(2*n/48)) for n in map(int, sys.argv[2].split(','))]
    else:
        designs = [tuple(map(int, s.split(','))) for s in sys.argv[1:]]
    for d in designs:
        lam, vp, gap, rel, tb, tt = measure(*d)
        print(f"(n,b,q,g,h)={d}: lambda={lam:.5f} (lambda - log2/2 = {lam-math.log(2)/2:.5f})  varpi={vp:.5f}  gap={gap:.5f}  "
              f"prec rel {rel:.1e}  [build {tb:.1f}s, total {tt:.1f}s, {time.time()-T0:.0f}s]", flush=True)
