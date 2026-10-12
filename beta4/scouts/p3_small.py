"""P3 data: for small primes, v_p of the determinant denominator vs the row bound (sum over rows of max v_p of the
entry denominators) and the largest single-entry v_p; B7 (beta4 #7) and C7 (Catalan, same rows). usage: python3 p3_small.py n"""
import sys, math
from fractions import Fraction as Fr
from nu_arrays import c
import scout_b4n as N
from scout_b4l import det_exact
n = int(sys.argv[1])
def vp(x, p):
    v = 0
    while x % p == 0: x //= p; v += 1
    return v
Mnu, MA, ML = N.arrays(2*n+6)
C7 = [[MA(2*r, 2*k).get('1', Fr(0)) + Fr(1, 3)*MA(2*r, 2*k).get('G', Fr(0)) for k in range(n)] for r in range(n)]
J = [2*k for k in range(n)]; B7 = [[None]*n for _ in range(n)]
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
    for k in range(n): B7[r][k] = col[2*k].get('1', Fr(0)) + Fr(1, 3)*col[2*k].get('b4', Fr(0))
for nm, M in (('C7', C7), ('B7', B7)):
    D = det_exact(M).denominator
    out = []
    for p in (2, 3, 5, 7, 11, 13):
        if p > 2*n: break
        rowb = sum(max(vp(x.denominator, p) for x in row) for row in M)
        mx = max(vp(x.denominator, p) for row in M for x in row)
        out.append(f"p={p}: det {vp(D, p)}, row bound {rowb}, max entry {mx}, log_p(2n) {math.log(2*n, p):.2f}")
    print(nm, f"n={n}:", " | ".join(out), flush=True)
