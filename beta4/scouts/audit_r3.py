"""Audit of route 3 + per-prime split of the determinant denominators.
A1  bases used by the even design, exact value vs quadrature at 100 digits: K-_0 = -8 beta4 (nu), K+_2 = -2 X3 (nu), KL+_2 = 4 - 4G (L), KL-_0 = X3 (L).
A2  whole entries of the overlap determinant (n = 4): value from the exact decomposition (rational + beta4 coefficient) vs
    an independent evaluation: rows/columns/auxiliary coefficients applied to quadrature values of the three arrays (30 digits).
D   per-prime split of -log den(Delta at 1/3)/n^2 for: beta4 overlap design, Catalan same design, half-integer weight 4 and weight 2."""
import math, mpmath as mp, sys
from fractions import Fraction as Fr
from math import comb
from nu_arrays import c
import scout_b4n as N
from scout_b4l import det_exact
mp.mp.dps = 100
pi = mp.pi; G = +mp.catalan; L2 = mp.log(2)
X3 = 8*mp.im(mp.polylog(3, mp.mpc(1, 1)/2)) - pi**3/16 - pi/4*L2**2
B4 = (mp.zeta(4, mp.mpf(1)/4) - mp.zeta(4, mp.mpf(3)/4))/256
def lm(th): return -mp.log(2)-2*mp.log(mp.sin(pi/4-th/2))
def lp(th): return -mp.log(1+mp.sin(th))
def ker(j, t, Lg, p):
    if abs(t) < mp.mpf('0.3'): return mp.nsum(lambda u: t**u/(u+j+1)**p, [0, mp.inf])
    head = Lg if p == 1 else mp.polylog(2, t)
    return (head - sum(t**k/mp.mpf(k)**p for k in range(1, j+1)))/t**(j+1)
def quadv(w, i, j, p):
    f = lambda th: w(th)*(mp.sin(th)**i*ker(j, mp.sin(th), lm(th), p) + (-mp.sin(th))**i*ker(j, -mp.sin(th), lp(th), p))
    return mp.re(mp.quad(f, [0, pi/8, pi/4, 3*pi/8, pi/2]))
wnu = lambda th: mp.sin(th)*(th*th - pi*th); wA = lambda th: mp.sin(th)
if 'A1' in sys.argv:
    for nm, v, ex in (('K-_0 (nu)', quadv(wnu, 0, 0, 1), -8*B4), ('K+_2 (nu)', quadv(wnu, 0, 2, 1), -2*X3),
                      ('KL+_2 (L)', quadv(wA, 0, 2, 2), 4-4*G), ('KL-_0 (L)', quadv(wA, 0, 0, 2), X3)):
        print('A1', nm, 'diff =', mp.nstr(v-ex, 3), flush=True)
if 'A2' in sys.argv:
    mp.mp.dps = 30
    n = 4; Mnu, MA, ML = N.arrays(20)
    cols = [{2*k: Fr(1)} for k in range(n)]; J = sorted(j for cc in cols for j in cc)
    qcache = {}
    def Q(nm, i, j):
        if (nm, i, j) not in qcache:
            qcache[(nm, i, j)] = quadv(wnu, i, j, 1) if nm == 'nu' else quadv(wA, i, j, 1 if nm == 'A' else 2)
        return qcache[(nm, i, j)]
    val = {'1': 1, 'G': G, 'X3': X3, 'b4': B4}
    for r in (0, 3):
        a = {2*r: Fr(1)}
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
        for k in range(n):
            j = 2*k; e = col[j]
            exact = sum(mp.mpf(v.numerator)/v.denominator*val[kk] for kk, v in e.items())
            indep = sum(mp.mpf(av.numerator)/av.denominator*Q('nu', i, j) for i, av in a.items()) \
                  + sum(mp.mpf(av.numerator)/av.denominator*Q('L', i2, j) for i2, av in aL.items()) \
                  + sum(mp.mpf(av.numerator)/av.denominator*Q('A', i2, j) for i2, av in aA.items())
            print('A2', (r, k), 'keys', sorted(e), 'diff =', mp.nstr(exact-indep, 3), flush=True)
def primes(N_):
    s = [1]*(N_+1); s[0] = s[1] = 0
    for i in range(2, int(N_**0.5)+1):
        if s[i]:
            for j in range(i*i, N_+1, i): s[j] = 0
    return [i for i in range(N_+1) if s[i]]
def split(den, n, P):
    out = {}; rest = den
    for p in P:
        v = 0
        while rest % p == 0: rest //= p; v += 1
        if v: out[p] = v
    return out, rest
if 'D' in sys.argv:
    P = primes(400)
    for n in (10, 12):
        Mnu, MA, ML = N.arrays(2*n+6)
        # Catalan same design
        E = [[MA(2*r, 2*k) for k in range(n)] for r in range(n)]
        Dc = det_exact([[e.get('1', Fr(0))+Fr(1, 3)*e.get('G', Fr(0)) for e in row] for row in E])
        # beta4 overlap design: rebuild via run() internals
        import io, contextlib
        # reuse N.run but we need the denominator: recompute matrix
        cols = [{2*k: Fr(1)} for k in range(n)]; J = sorted(j for cc in cols for j in cc)
        Frat = [[None]*n for _ in range(n)]
        for r in range(n):
            a = {2*r: Fr(1)}; col = {j: {} for j in J}
            for j in J: col[j] = N.addd(col[j], Mnu(2*r, j))
            aL = {}
            for j in sorted(J, reverse=True): aL[j] = -col[j].get('X3', Fr(0)) - sum((aL[i2]*c((i2-j)//2) for i2 in aL if i2 > j), Fr(0))
            for j in J:
                for i2, av in aL.items(): col[j] = N.addd(col[j], ML(i2, j), av)
            aA = {}
            for j in sorted(J, reverse=True): aA[j] = (-col[j].get('G', Fr(0)) - sum((aA[i2]*4*c((i2-j)//2) for i2 in aA if i2 > j), Fr(0)))/4
            for j in J:
                for i2, av in aA.items(): col[j] = N.addd(col[j], MA(i2, j), av)
            for k in range(n): Frat[r][k] = col[2*k].get('1', Fr(0)) + Fr(1, 3)*col[2*k].get('b4', Fr(0))
        Db = det_exact(Frat)
        for nm, D_ in (('Catalan same design', Dc), ('beta4 overlap', Db)):
            sp, rest = split(D_.denominator, n, P)
            groups = {'p=2': 0.0, 'p=3': 0.0, '5<=p<=n': 0.0, 'n<p<=2n': 0.0, '2n<p<=4n': 0.0, 'p>4n': 0.0}
            for p, v in sp.items():
                g = 'p=2' if p == 2 else 'p=3' if p == 3 else '5<=p<=n' if p <= n else 'n<p<=2n' if p <= 2*n else '2n<p<=4n' if p <= 4*n else 'p>4n'
                groups[g] += v*math.log(p)/n**2
            print(f"D n={n} {nm:20s}: total {-math.log(D_.denominator)/n**2:.3f}; by primes (per n^2): " + ", ".join(f"{g} {-x:.3f}" for g, x in groups.items()) + (f"; unfactored {rest}" if rest > 1 else ""), flush=True)
