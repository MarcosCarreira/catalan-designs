"""Hankel determinants of Rivoal-Zudilin moments (scout).
mu_j(N) = sum_{t>=0} (-1)^t R_N(t) u^{2j}, u = t + h0/2, R_N the very-well-poised q=6 rational function of rz_h1.py (parameter eta).
The factor u^{2j} keeps the symmetry t -> -h0-t, so mu_j = c0 + c2 beta(2) + c4 beta(4) (odd beta's vanish; checked).
H_m = det(mu_{r+k})_{r,k<m}. Measures lambda = log|H_m|/m^2 and varpi = -log den(H_m with beta(2), beta(4) -> 1/3)/m^2.
usage: python3 rz_hankel.py eta(comma) N m"""
import sys, math, time
from fractions import Fraction as F
import mpmath as mp
exec(open('rz_scout.py').read().split("etas=")[0])
T0 = time.time()
def moment(eta, n, jpow):
    q = len(eta)-1; e0 = eta[0]; h0 = e0*n+1; hq = eta[q]*n+1
    h = [None]+[F(2*eta[j]*n+1, 2) for j in range(1, q)]+[hq]
    gam = F(1)
    for j in range(1, q): gam *= math.factorial((e0-2*eta[j])*n)
    gam /= math.factorial(eta[q]*n)**2
    num = [F(h0)*gam, 2*gam]
    for i in range(hq-1): num = pmul(num, [F(1+i), F(1)])
    for i in range(hq-1): num = pmul(num, [F(1+h0-hq+i), F(1)])
    for _ in range(2*jpow): num = pmul(num, [F(h0, 2), F(1)])
    poles = {}
    for j in range(1, q):
        L = 1+h0-int(2*h[j])
        for i in range(L):
            p = -(h[j]+i); poles[p] = poles.get(p, 0)+1
    roots = []
    for p, ml in poles.items(): roots += [p]*ml
    assert len(num)-1 < len(roots), "degree too high"
    mm = max(poles.values()); coef = {}
    for p, ml in poles.items():
        den = [F(1)]
        for r in roots:
            if r != p: den = pmul(den, [-r, F(1)])
        nu = [F(0)]*(mm-ml)+shift(num, p); de = shift(den, p)
        ser = series_div(nu, de, mm); k = -p+F(1, 2)
        for j in range(1, mm+1): coef[(k, j)] = ser[mm-j]
    cb = {j: F(0) for j in range(1, mm+1)}; c0 = F(0)
    for (k, j), c in coef.items():
        k = int(k); sg = (-1)**(k+1)
        cb[j] += c*sg*2**j
        c0 -= c*(-1)**k*sum(F((-1)**l)/F(2*l-1, 2)**j for l in range(1, k))
    assert all(cb[j] == 0 for j in cb if j % 2 == 1), cb
    assert all(cb[j] == 0 for j in cb if j > 4), cb
    return c0, cb.get(2, F(0)), cb.get(4, F(0)), len(roots)-(len(num)-1)
def det_exact(M):
    n = len(M); scale = 1; A = []
    for row in M:
        Lc = 1
        for x in row: Lc = Lc*x.denominator//math.gcd(Lc, x.denominator)
        scale *= Lc; A.append([int(x*Lc) for x in row])
    sign = 1; prev = 1
    for cc in range(n-1):
        p = next((r for r in range(cc, n) if A[r][cc] != 0), None)
        if p is None: return F(0)
        if p != cc: A[cc], A[p] = A[p], A[cc]; sign = -sign
        for r in range(cc+1, n):
            for k in range(cc+1, n): A[r][k] = (A[r][k]*A[cc][cc] - A[r][cc]*A[cc][k])//prev
            A[r][cc] = 0
        prev = A[cc][cc]
    return F(sign*A[n-1][n-1], scale)
if __name__ == "__main__":
    eta = tuple(int(x) for x in sys.argv[1].split(',')); N = int(sys.argv[2]); m = int(sys.argv[3])
    mus = [moment(eta, N, j) for j in range(2*m-1)]
    print(f"eta={eta} N={N} m={m}: degree excess at j=2m-2: {mus[-1][3]}  [{time.time()-T0:.0f}s]", flush=True)
    mp.mp.dps = 80 + 3*m*m + 40*N
    b2 = +mp.catalan; b4 = (mp.zeta(4, mp.mpf(1)/4)-mp.zeta(4, mp.mpf(3)/4))/256
    fv = lambda x: mp.mpf(x.numerator)/x.denominator
    vals = [fv(c0)+fv(c2)*b2+fv(c4)*b4 for c0, c2, c4, _ in mus]
    if 'check' in sys.argv:   # positive control: mu_j by direct summation
        q = len(eta)-1; e0 = eta[0]; h0 = e0*N+1
        hh = [h0]+[mp.mpf(2*eta[j]*N+1)/2 for j in range(1, q)]+[eta[q]*N+1]
        g = mp.mpf(1)
        for j in range(1, q): g *= mp.factorial((e0-2*eta[j])*N)
        g /= mp.factorial(eta[q]*N)**2
        for jj in (0, 2*m-2):
            def term(t):
                v = (2*t+h0)*mp.gamma(t+h0)/mp.gamma(t+1)*(t+mp.mpf(h0)/2)**(2*jj)
                for j in range(1, q+1): v *= mp.gamma(t+hh[j])/mp.gamma(t+1+h0-hh[j])
                return g*v*(-1)**int(t)
            print(f"  control j={jj}: exact {mp.nstr(vals[jj], 20)} direct {mp.nstr(mp.nsum(term, [0, mp.inf]), 20)}", flush=True)
    H = mp.matrix(m, m)
    for r in range(m):
        for k in range(m): H[r, k] = vals[r+k]
    Ht = mp.det(H)
    Dh = det_exact([[mus[r+k][0]+F(1, 3)*(mus[r+k][1]+mus[r+k][2]) for k in range(m)] for r in range(m)])
    lam = float(mp.log(abs(Ht)))/m**2; vp = -math.log(Dh.denominator)/m**2
    sizes = [float(mp.log(abs(v))) for v in vals]
    print(f"  single moments: log|mu_0| = {sizes[0]:.2f}, log|mu_{2*m-2}| = {sizes[-1]:.2f}; log den(mu_0) = {math.log(math.lcm(*(x.denominator for x in mus[0][:3]))):.2f}")
    print(f"  Hankel: lambda = {lam:.4f}  varpi = {vp:.4f}  gap = {lam-vp:.4f}   (log|H| = {lam*m*m:.1f}, log den = {-vp*m*m:.1f})  [{time.time()-T0:.0f}s]", flush=True)
