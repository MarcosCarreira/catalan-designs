"""Content-normalized measure for the arcsine-family determinants: PRE (OpenAI, ratios 48:7:4:4:2), C7 (Catalan, arcsine
measure, even indices), B7 (beta(4), weighted measure #7). Q = primitive integer polynomial of det(Frat + X Fc); value at C."""
import sys, math, time
from fractions import Fraction as Fr
from flint import fmpq, fmpq_poly, fmpz_mat
import mpmath as mp
sys.path.insert(0, '.')
from nu_arrays import c
def mats(name, n):
    if name == 'PRE':
        import full_cat as FC
        d = (n, round(7*n/48), round(4*n/48), round(4*n/48), round(2*n/48))
        return FC.build(*d)
    import scout_b4n as N
    Mnu, MA, ML = N.arrays(2*n+6)
    if name == 'C7':
        E = [[MA(2*r, 2*k) for k in range(n)] for r in range(n)]
        return [[e.get('1', Fr(0)) for e in row] for row in E], [[e.get('G', Fr(0)) for e in row] for row in E]
    J = [2*k for k in range(n)]; A = [[None]*n for _ in range(n)]; B = [[None]*n for _ in range(n)]
    for r in range(n):
        col = {j: N.addd({}, Mnu(2*r, j)) for j in J}
        aL = {}
        for j in sorted(J, reverse=True): aL[j] = -col[j].get('X3', Fr(0)) - sum((aL[i2]*c((i2-j)//2) for i2 in aL if i2 > j), Fr(0))
        for j in J:
            for i2, av in aL.items(): col[j] = N.addd(col[j], ML(i2, j), av)
        aA = {}
        for j in sorted(J, reverse=True): aA[j] = (-col[j].get('G', Fr(0)) - sum((aA[i2]*4*c((i2-j)//2) for i2 in aA if i2 > j), Fr(0)))/4
        for j in J:
            for i2, av in aA.items(): col[j] = N.addd(col[j], MA(i2, j), av)
        for k in range(n):
            assert all(kk in ('1', 'b4') for kk in col[2*k]); A[r][k] = col[2*k].get('1', Fr(0)); B[r][k] = col[2*k].get('b4', Fr(0))
    return A, B
def measure(name, n):
    T0 = time.time()
    A, B = mats(name, n)
    L = 1
    for row in A + B:
        for v in row: L = L*v.denominator//math.gcd(L, v.denominator)
    Ai = [[int(v*L) for v in row] for row in A]; Bi = [[int(v*L) for v in row] for row in B]
    xs = list(range(n+1)); ys = []
    for x in xs:
        ys.append(fmpq(int(fmpz_mat(n, n, [Ai[r][k] + x*Bi[r][k] for r in range(n) for k in range(n)]).det())))
    coef = ys[:]
    for k in range(1, n+1):
        for i in range(n, k-1, -1): coef[i] = (coef[i]-coef[i-1])/(xs[i]-xs[i-k])
    t = fmpq_poly([0, 1]); P = fmpq_poly([coef[n]])
    for k in range(n-1, -1, -1): P = P*(t-xs[k]) + coef[k]
    deg = P.degree()
    g = 0; D = 1
    for k in range(deg+1): g = math.gcd(g, int(P[k].p)); D = D*int(P[k].q)//math.gcd(D, int(P[k].q))
    Q = P*fmpq(D, g)
    bits = max(abs(int(Q[k].p)).bit_length() for k in range(deg+1))
    d = int(1.3*0.302*bits) + 300
    for _ in range(5):
        vals = []
        for dd in (d, d+1000):
            mp.mp.dps = dd
            C = +mp.catalan if name in ('PRE', 'C7') else (mp.zeta(4, mp.mpf(1)/4)-mp.zeta(4, mp.mpf(3)/4))/256
            v = mp.mpf(0)
            for k in range(deg, -1, -1): v = v*C + int(Q[k].p)
            vals.append(v)
        if abs(vals[0]-vals[1]) <= abs(vals[1])*mp.mpf(10)**-30: break
        d *= 2
    lq = float(mp.log(abs(vals[1])))
    print(f"{name} n={n} deg={deg}: log|Q(C)| = {lq:9.2f}  per n^2 = {lq/n**2:+.4f}  log b-bound = {-lq/deg:+.3f}  [{time.time()-T0:.1f}s]", flush=True)
if __name__ == "__main__":
    for name in sys.argv[1].split(','):
        for n in map(int, sys.argv[2].split(',')): measure(name, n)
