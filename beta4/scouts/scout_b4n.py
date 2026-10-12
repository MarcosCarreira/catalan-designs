"""Route 3, overlap design (scout). Arrays on even indices:
  nu-array  M(i,j)  = int dnu t^i int s^j/(1-ts)            constants beta4, G, X3   (nu_arrays.py + K+ even recurrence)
  A-array   MA(i,j) = int dmu t^i int s^j/(1-ts)            constants G              (mu = |t| dt/sqrt(1-t^2))
  L-array   L(i,j)  = int dmu t^i int s^j(-log s)/(1-ts)    constants X3, G
Rows t^{2(alpha+r)}(1-t^2)^h, columns s^{2(b+k)}(1-s^2)^q. Auxiliary L-rows on the column support J cancel X3,
then auxiliary A-rows on J cancel G (both triangular). Entries in Q + Q beta4.
Recurrences (route-3 notes): KL-_d = (1-1/d) KL-_{d-2} + (KA-_{d-2} - kp_2 + sum_{v<d-2} m_v/(v+3))/d^2, KL-_0 = X3;
KL+_{d+2} = d/(d+1) [KL+_d - m_0/(d+1)^2 + (kp_2 - kp_{d+2})/d^2], KL+_2 = 4 - 4G (PSLQ);
K+_{j+2} = [j K+_j + K+_2 - j F_j]/(j+1), F_j = -2[X3/j - 4G/j^2 + kp_j/j^2], K+_2 = -2 X3, K+_0 = -8 beta4.
usage: python3 scout_b4n.py n1,n2,.. rh rq ralpha rb  [check]"""
import sys, math, time
from fractions import Fraction as Fr
from math import comb
import mpmath as mp
from nu_arrays import build_nu, M_nu, c
from scout_b4l import det_exact, kp_arr
T0 = time.time()
def addd(a, b, f=Fr(1)):
    o = dict(a)
    for k, v in b.items(): o[k] = o.get(k, Fr(0)) + f*v
    return {k: v for k, v in o.items() if v != 0}
def arrays(top):
    m, km, nu, K = build_nu(top+2); kp = kp_arr(top+2)
    S3 = [Fr(0)]*(top+4)
    for d in range(3, top+4): S3[d] = S3[d-1] + m[d-3]/d
    KAm = lambda d: addd({'1': km[d]}, {'G': 4*c(d//2)} if d % 2 == 0 else {})
    KLm = {0: {'X3': Fr(1)}}
    for d in range(2, top+1, 2):
        KLm[d] = addd({k: (1-Fr(1, d))*v for k, v in KLm[d-2].items()}, addd(KAm(d-2), {'1': -kp[2] + S3[d]}), Fr(1, d*d))
    KLp = {2: {'1': Fr(4), 'G': Fr(-4)}}
    for d in range(2, top-1, 2):
        KLp[d+2] = {k: Fr(d, d+1)*v for k, v in addd(KLp[d], {'1': -m[0]/(d+1)**2 + (kp[2]-kp[d+2])/d**2}).items()}
    Kp = {0: {'b4': Fr(-8)}, 2: {'X3': Fr(-2)}}
    for j in range(2, top-1, 2):
        F = {'X3': Fr(-2, j), 'G': Fr(8, j*j), '1': -2*kp[j]/j**2}
        tot = addd({k: j*v for k, v in Kp[j].items()}, Kp[2]); tot = addd(tot, F, Fr(-j))
        Kp[j+2] = {k: v/(j+1) for k, v in tot.items()}
    def Mnu(i, j):
        if i >= j: return M_nu(i, j, nu, K)
        return addd(Kp[j-i], {'1': -sum((nu[l-1]/(j-i+l) for l in range(1, i+1)), Fr(0))})
    def MA(i, j):
        if i >= j: return addd(KAm(i-j), {'1': -sum((m[i-j+l-1]/l for l in range(1, j+1)), Fr(0))})
        return {'1': kp[j-i] - sum((m[l-1]/(j-i+l) for l in range(1, i+1)), Fr(0))}
    def ML(i, j):
        if i >= j: return addd(KLm[i-j], {'1': -sum((m[i-j+l-1]/Fr(l*l) for l in range(1, j+1)), Fr(0))})
        return addd(KLp[j-i], {'1': -sum((m[l-1]/Fr((j-i+l)**2) for l in range(1, i+1)), Fr(0))})
    return Mnu, MA, ML
def run(n, h, q, alpha, b, check=False):
    top = 2*max(alpha+n-1+h, b+n-1+q) + 4
    Mnu, MA, ML = arrays(top)
    cols = [{2*(b+k+l): Fr((-1)**l*comb(q, l)) for l in range(q+1)} for k in range(n)]
    J = sorted(set(j for cc in cols for j in cc))
    Frat = [[Fr(0)]*n for _ in range(n)]; Fb4 = [[Fr(0)]*n for _ in range(n)]; resid = 0
    for r in range(n):
        a = {2*(alpha+r+l): Fr((-1)**l*comb(h, l)) for l in range(h+1)}
        col = {j: {} for j in J}
        for j in J:
            for i, av in a.items(): col[j] = addd(col[j], Mnu(i, j), av)
        aL = {}
        for j in sorted(J, reverse=True):
            s = -col[j].get('X3', Fr(0)) - sum((aL[i2]*c((i2-j)//2) for i2 in aL if i2 > j), Fr(0))
            aL[j] = s
        for j in J:
            for i2, av in aL.items(): col[j] = addd(col[j], ML(i2, j), av)
        aA = {}
        for j in sorted(J, reverse=True):
            s = -col[j].get('G', Fr(0)) - sum((aA[i2]*4*c((i2-j)//2) for i2 in aA if i2 > j), Fr(0))
            aA[j] = s/4
        for j in J:
            for i2, av in aA.items(): col[j] = addd(col[j], MA(i2, j), av)
        for k, ck in enumerate(cols):
            e = {}
            for j, cv in ck.items(): e = addd(e, col[j], cv)
            if any(k2 not in ('1', 'b4') for k2 in e): resid += 1
            Frat[r][k] = e.get('1', Fr(0)); Fb4[r][k] = e.get('b4', Fr(0))
    Dh = det_exact([[Frat[r][k] + Fr(1, 3)*Fb4[r][k] for k in range(n)] for r in range(n)])
    outs = []
    for dps in (60 + 6*n*n, 120 + 6*n*n):
        mp.mp.dps = dps
        B4 = (mp.zeta(4, mp.mpf(1)/4) - mp.zeta(4, mp.mpf(3)/4))/256
        Mm = mp.matrix(n, n)
        for r in range(n):
            for k in range(n):
                Mm[r, k] = mp.mpf(Frat[r][k].numerator)/Frat[r][k].denominator + B4*mp.mpf(Fb4[r][k].numerator)/Fb4[r][k].denominator
        outs.append(mp.det(Mm))
    rel = abs(outs[0]-outs[1])/abs(outs[1])
    lam = float(mp.log(abs(outs[1]))/n**2); vp = -math.log(Dh.denominator)/n**2
    return resid, lam, vp, lam-vp, float(rel), Dh != 0
if __name__ == "__main__":
    if len(sys.argv) > 6 and sys.argv[6] == 'check':
        # positive control: exact L, A, nu entries against quadrature (40 digits)
        mp.mp.dps = 40
        pi = mp.pi; G = +mp.catalan; X3 = 8*mp.im(mp.polylog(3, mp.mpc(1, 1)/2)) - pi**3/16 - pi/4*mp.log(2)**2
        B4 = (mp.zeta(4, mp.mpf(1)/4) - mp.zeta(4, mp.mpf(3)/4))/256
        val = {'1': 1, 'G': G, 'X3': X3, 'b4': B4}
        Mnu, MA, ML = arrays(30)
        def lm(th): return -mp.log(2)-2*mp.log(mp.sin(pi/4-th/2))
        def lp(th): return -mp.log(1+mp.sin(th))
        def ker(j, t, Lg, p):   # sum_u t^u/(u+j+1)^p
            if abs(t) < mp.mpf('0.3'): return mp.nsum(lambda u: t**u/(u+j+1)**p, [0, mp.inf])
            head = Lg if p == 1 else mp.polylog(2, t)
            return (head - sum(t**k/mp.mpf(k)**p for k in range(1, j+1)))/t**(j+1)
        def quadv(w, i, j, p):
            f = lambda th: w(th)*(mp.sin(th)**i*ker(j, mp.sin(th), lm(th), p) + (-mp.sin(th))**i*ker(j, -mp.sin(th), lp(th), p))
            return mp.re(mp.quad(f, [0, pi/4, pi/2]))
        wnu = lambda th: mp.sin(th)*(th*th - pi*th); wA = lambda th: mp.sin(th)
        worst = 0
        for (i, j) in [(0, 2), (2, 6), (4, 8), (6, 2), (8, 8), (2, 10)]:
            for nm, F_, w, p in (('nu', Mnu, wnu, 1), ('A', MA, wA, 1), ('L', ML, wA, 2)):
                ex = sum(mp.mpf(v.numerator)/v.denominator*val[k] for k, v in F_(i, j).items())
                d = abs(ex - quadv(w, i, j, p)); worst = max(worst, d)
                print(nm, (i, j), mp.nstr(d, 3), {k: str(v) for k, v in F_(i, j).items() if k != '1'}, flush=True)
        print('worst', mp.nstr(worst, 3)); sys.exit()
    ns = [int(x) for x in sys.argv[1].split(',')]; rh, rq, ra, rb = (float(x) for x in sys.argv[2:6])
    for n in ns:
        h, q, al, b = round(rh*n), round(rq*n), round(ra*n), round(rb*n)
        resid, lam, vp, gap, rel, nz = run(n, h, q, al, b)
        print(f"n={n} (h,q,alpha,b)=({h},{q},{al},{b}): residues {resid}; lambda={lam:.4f} varpi={vp:.4f} gap={gap:.4f} prec {rel:.1e} [{time.time()-T0:.0f}s]", flush=True)
