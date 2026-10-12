"""Route 3 determinant (scout). Measure nu (nu_arrays.py); rows t^{2(alpha+r)} (1-t^2)^h with alpha = b+n-1+q+da,
so every row exponent >= every column index; columns s^{2(b+k)} (1-s^2)^q; extra arcsine-array rows a'_r on the
column support J cancel G exactly (triangular system). Entries in Q + Q beta(4).
Measures lambda = log|Delta|/n^2, varpi = -log den(Delta with beta4 -> 1/3)/n^2, gap = lambda - varpi.
usage: python3 scout_b4l.py n1,n2,... h q da b"""
import sys, math, time
from fractions import Fraction as Fr
from math import comb
import mpmath as mp
from nu_arrays import build_nu, M_nu, c
T0 = time.time()
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
def kp_arr(top):
    kp = [Fr(0)]*(top+3)
    for d in range(2, top+3): kp[d] = ((d-2)*kp[d-2] + Fr(2, d-1))/(d-1)
    return kp
def run(n, h, q, da, b):
    alpha = b + n - 1 + q + da
    I = 2*(alpha + n - 1 + h)
    m, km, nu, K = build_nu(I + 2)
    kp = kp_arr(I + 2)
    cols = [{2*(b+k+l): Fr((-1)**l*comb(q, l)) for l in range(q+1)} for k in range(n)]
    J = sorted(set(j for cc in cols for j in cc))
    def MA(i, j):   # arcsine array, even i, j: (rational part, G coefficient)
        if i >= j:
            d = i-j; rat = km[d] - sum((m[d+l-1]/l for l in range(1, j+1)), Fr(0))
            return rat, 4*c(d//2)
        d = j-i; rat = kp[d] - sum((m[l-1]/(d+l) for l in range(1, i+1)), Fr(0))
        return rat, Fr(0)
    Frat = [[Fr(0)]*n for _ in range(n)]; Fb4 = [[Fr(0)]*n for _ in range(n)]; resid = 0
    for r in range(n):
        a = {2*(alpha+r+l): Fr((-1)**l*comb(h, l)) for l in range(h+1)}
        # nu part per column index j
        nuj = {}
        for j in J:
            tot = {}
            for i, av in a.items():
                for key, v in M_nu(i, j, nu, K).items(): tot[key] = tot.get(key, Fr(0)) + av*v
            nuj[j] = tot
        # auxiliary arcsine rows: sum_{i' in J, i' >= j} a'_{i'} 4 c((i'-j)/2) = -G_nu(j)
        ap = {}
        for j in sorted(J, reverse=True):
            s = -nuj[j].get('G', Fr(0)) - sum((ap[i2]*4*c((i2-j)//2) for i2 in ap if i2 > j), Fr(0))
            ap[j] = s/4
        for k, ck in enumerate(cols):
            rat = Fr(0); g = Fr(0); b4 = Fr(0)
            for j, cv in ck.items():
                rat += cv*nuj[j].get('1', Fr(0)); g += cv*nuj[j].get('G', Fr(0)); b4 += cv*nuj[j].get('b4', Fr(0))
                for i2, av in ap.items():
                    ra, ga = MA(i2, j); rat += cv*av*ra; g += cv*av*ga
            if g != 0: resid += 1
            Frat[r][k] = rat; Fb4[r][k] = b4
    Dh = det_exact([[Frat[r][k] + Fr(1, 3)*Fb4[r][k] for k in range(n)] for r in range(n)])
    outs = []
    for dps in (60 + 6*n*n, 120 + 6*n*n):
        mp.mp.dps = dps
        B4 = (mp.zeta(4, mp.mpf(1)/4) - mp.zeta(4, mp.mpf(3)/4))/256
        Mm = mp.matrix(n, n)
        for r in range(n):
            for k in range(n):
                x, y = Frat[r][k], Fb4[r][k]
                Mm[r, k] = mp.mpf(x.numerator)/x.denominator + B4*mp.mpf(y.numerator)/y.denominator
        outs.append(mp.det(Mm))
    rel = abs(outs[0]-outs[1])/abs(outs[1]) if outs[1] != 0 else mp.inf
    lam = float(mp.log(abs(outs[1]))/n**2) if outs[1] != 0 else None
    vp = -math.log(Dh.denominator)/n**2 if Dh != 0 else None
    return dict(alpha=alpha, resid=resid, lam=lam, varpi=vp, gap=(lam-vp) if (lam is not None and vp is not None) else None, rel=float(rel))
if __name__ == "__main__":
    ns = [int(x) for x in sys.argv[1].split(',')]; rh, rq, rda, rb = (float(x) for x in sys.argv[2:6])
    for n in ns:
        h, q, da, b = round(rh*n), round(rq*n), round(rda*n), round(rb*n)
        R = run(n, h, q, da, b)
        print(f"n={n} (h,q,da,b)=({h},{q},{da},{b}) alpha={R['alpha']}: G residues {R['resid']}; lambda={R['lam']:.4f} varpi={R['varpi']:.4f} "
              f"gap={R['gap']:.4f}  prec rel {R['rel']:.1e}  [{time.time()-T0:.0f}s]", flush=True)
