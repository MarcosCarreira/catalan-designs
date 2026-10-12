"""P1: per-prime profile of the exact determinant denominators (C -> 1/3), binned by p/n.
Constructions: B7 = beta(4) overlap design (#7), C7 = Catalan same rows/columns (#7c), PRE = preprint's determinant (full_cat.py),
H4/H2 = half-integer family weight 4 / weight 2 (b4_balance.py, design (0,0,0,0)).
usage: python3 p1_profile.py <name> n1,n2,..   (name in B7 C7 PRE H4 H2)"""
import sys, math, time
from fractions import Fraction as Fr
from nu_arrays import c
T0 = time.time()
def primes(N_):
    s = [1]*(N_+1); s[0] = s[1] = 0
    for i in range(2, int(N_**0.5)+1):
        if s[i]:
            for j in range(i*i, N_+1, i): s[j] = 0
    return [i for i in range(N_+1) if s[i]]
def den_B7_C7(n, which):
    import scout_b4n as N
    from scout_b4l import det_exact
    Mnu, MA, ML = N.arrays(2*n+6)
    if which == 'C7':
        E = [[MA(2*r, 2*k) for k in range(n)] for r in range(n)]
        return det_exact([[e.get('1', Fr(0))+Fr(1, 3)*e.get('G', Fr(0)) for e in row] for row in E]).denominator
    J = [2*k for k in range(n)]; Fm = [[None]*n for _ in range(n)]
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
        for k in range(n): Fm[r][k] = col[2*k].get('1', Fr(0)) + Fr(1, 3)*col[2*k].get('b4', Fr(0))
    return det_exact(Fm).denominator
def den_PRE(n):
    import full_cat as FC
    d = (n, round(7*n/48), round(4*n/48), round(4*n/48), round(2*n/48))
    Frat, FG = FC.build(*d)
    return FC.det_exact([[Frat[r][k] + Fr(1, 3)*FG[r][k] for k in range(n)] for r in range(n)]).denominator
def den_H(n, W):
    sys.argv = ['x']
    src = open('b4_balance.py').read().split('if __name__')[0]
    g = {'__file__': 'b4_balance.py', '__name__': 'b4b'}; exec(src, g)
    Frat, Ft, resid = g['build'](n, 0, 0, 0, 0, W)
    return g['det_exact']([[Frat[r][k] + Fr(1, 3)*Ft[r][k] for k in range(n)] for r in range(n)]).denominator
name = sys.argv[1]
for n in map(int, sys.argv[2].split(',')):
    den = {'B7': lambda: den_B7_C7(n, 'B7'), 'C7': lambda: den_B7_C7(n, 'C7'), 'PRE': lambda: den_PRE(n),
           'H4': lambda: den_H(n, 4), 'H2': lambda: den_H(n, 2)}[name]()
    P = primes(8*n+20); bins = {}; rest = den; maxp = 0
    for p in P:
        v = 0
        while rest % p == 0: rest //= p; v += 1
        if v:
            b = min(int(4*p/n), 15)  # bins of width 0.25 in p/n
            bins[b] = bins.get(b, 0.0) + v*math.log(p)/n**2; maxp = p
    tot = math.log(den)/n**2
    line = " ".join(f"[{b/4:.2f},{(b+1)/4:.2f}):{bins.get(b, 0):.2f}" for b in range(0, max(bins)+1))
    print(f"{name} n={n}: total {tot:.3f} (largest prime {maxp}, {maxp/n:.2f} n){' UNFACTORED '+str(rest) if rest > 1 else ''}\n   {line}  [{time.time()-T0:.0f}s]", flush=True)
