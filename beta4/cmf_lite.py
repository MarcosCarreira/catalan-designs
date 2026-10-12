"""cmf_lite.py -- pFq conservative matrix field without external packages (Fractions + Python integers + mpmath).
Same conventions as ramanujantools 0.0.1 (DFinite / pFq): state matrices from the theta-equation
  theta prod_i (theta + y_i - 1) - z prod_i (theta + x_i) = 0 (made monic), C = companion matrix (entries -c_i in the last column,
  ones on the subdiagonal); unit step in x_i: C/x_i + I; unit step in y_i: (C(y_i+1)/y_i + I)^(-1).
Trajectory step T(P) = product of unit steps from P to P + t; walk W_N = T(P0) T(P0+t) ... T(P0+(N-1)t);
p_N = W[0][-1], q_N = W[1][-1]; delta = -(1 + log|L - p/q| / log q) with p/q reduced and L from depth 2N."""
from fractions import Fraction as F
from math import gcd
import mpmath as mp
def pmul(a, b):
    o = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): o[i+j] += x*y
    return o
def companion(xs, ys, z):
    left = [F(0), F(1)]                                   # theta
    for y in ys: left = pmul(left, [y-1, F(1)])
    right = [F(1)]
    for x in xs: right = pmul(right, [x, F(1)])
    deg = max(len(left), len(right))
    left += [F(0)]*(deg-len(left)); right += [F(0)]*(deg-len(right))
    poly = [l - z*r for l, r in zip(left, right)]
    while poly and poly[-1] == 0: poly.pop()
    lc = poly[-1]; c = [v/lc for v in poly[:-1]]; s = len(c)
    C = [[F(0)]*s for _ in range(s)]
    for i in range(s):
        C[i][s-1] = -c[i]
        if i >= 1: C[i][i-1] = F(1)
    return C
def inv(A):
    n = len(A); M = [row[:] + [F(int(i == j)) for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0); M[c], M[p] = M[p], M[c]
        pv = M[c][c]; M[c] = [v/pv for v in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]; M[r] = [a - f*b for a, b in zip(M[r], M[c])]
    return [row[n:] for row in M]
def unit(xs, ys, z, kind, i):
    """unit step matrix at the point (xs, ys) along axis kind ('x' or 'y') number i, as an integer matrix up to a scalar"""
    if kind == 'x':
        C = companion(xs, ys, z); s = len(C)
        M = [[C[r][c]/xs[i] + (1 if r == c else 0) for c in range(s)] for r in range(s)]
    else:
        ys2 = list(ys); ys2[i] = ys[i] + 1
        C = companion(xs, ys2, z); s = len(C)
        M = inv([[C[r][c]/ys[i] + (1 if r == c else 0) for c in range(s)] for r in range(s)])
    L = 1
    for row in M:
        for v in row: L = L*v.denominator//gcd(L, v.denominator)
    return [[int(v*L) for v in row] for row in M]
def matmul(A, B):
    n, m, k = len(A), len(B[0]), len(B)
    return [[sum(A[i][l]*B[l][j] for l in range(k)) for j in range(m)] for i in range(n)]
def reduce(A):
    g = 0
    for row in A:
        for v in row: g = gcd(g, v)
    return [[v//g for v in row] for row in A] if g > 1 else A
def walk(xs0, ys0, z, tx, ty, depths):
    """returns {depth: (p, q)} for the requested depths (sorted)"""
    xs = [F(v) for v in xs0]; ys = [F(v) for v in ys0]; z = F(z)
    s = max(len(ys)+1, len(xs)); W = [[int(i == j) for j in range(s)] for i in range(s)]
    out = {}; depths = sorted(depths)
    import random
    rng = random.Random(12345)
    steps = [('y', i) for i in range(len(ty)) for _ in range(ty[i])] + [('x', i) for i in range(len(tx)) for _ in range(tx[i])]
    orders = [steps, steps[::-1]] + [rng.sample(steps, len(steps)) for _ in range(200)]
    for nstep in range(1, depths[-1]+1):
        for order in orders:            # the field is conservative: every monotone path gives the same product;
            try:                        # a path is skipped if it meets a singular unit matrix
                T = None; xs2 = list(xs); ys2 = list(ys)
                for kind, i in order:
                    for _ in range(1):
                        U = unit(xs2, ys2, z, kind, i)
                        T = U if T is None else matmul(T, U)
                        if kind == 'x': xs2[i] += 1
                        else: ys2[i] += 1
                break
            except (StopIteration, ZeroDivisionError):
                T = 'singular'
        if T == 'singular': raise ArithmeticError('all paths singular')
        W = reduce(matmul(W, T)); xs, ys = xs2, ys2
        if nstep in depths: out[nstep] = (W[0][-1], W[1][-1])
    return out
def delta(xs0, ys0, z, tx, ty, depth):
    r = walk(xs0, ys0, z, tx, ty, [depth, 2*depth])
    p, q = r[depth]; P, Q = r[2*depth]
    if q == 0 or Q == 0: return float('-inf')
    g = gcd(p, q); p //= g; q //= g
    if q < 0: p, q = -p, -q
    if abs(q) == 1: return float('inf')
    old = mp.mp.dps
    mp.mp.dps = int(3*len(str(abs(Q)))) + 50
    L = mp.mpf(P)/Q; err = abs(L - mp.mpf(p)/q)
    d = float(-(1 + mp.log(err)/mp.log(q))) if err != 0 else float('inf')
    mp.mp.dps = old
    return d
def apery_delta(depth):
    """Apery's recurrence for zeta(3): (n+1)^3 u_{n+1} = (34n^3+51n^2+27n+5) u_n - n^3 u_{n-1}; a: 0, 6; b: 1, 5."""
    a = [F(0), F(6)]; b = [F(1), F(5)]
    for n in range(1, 2*depth):
        c = 34*n**3+51*n**2+27*n+5
        a.append((c*a[n] - n**3*a[n-1])/(n+1)**3); b.append((c*b[n] - n**3*b[n-1])/(n+1)**3)
    r = a[depth]/b[depth]; p, q = r.numerator, r.denominator
    old = mp.mp.dps; mp.mp.dps = 3*len(str(q)) + 50
    L = mp.zeta(3); d = float(-(1 + mp.log(abs(L - mp.mpf(p)/q))/mp.log(q))); mp.mp.dps = old
    return d
def deltas(xs0, ys0, z, tx, ty, depths):
    """delta at each depth d in depths (limit taken from depth 2*max(depths)); one walk."""
    D = sorted(set(list(depths) + [2*max(depths)]))
    r = walk(xs0, ys0, z, tx, ty, D)
    P, Q = r[2*max(depths)]
    out = []
    for d in depths:
        p, q = r[d]
        if q == 0 or Q == 0: out.append(float('-inf')); continue
        g = gcd(p, q); p //= g; q //= g
        if q < 0: p, q = -p, -q
        if q == 1: out.append(float('inf')); continue
        old = mp.mp.dps; mp.mp.dps = int(3*len(str(abs(Q)))) + 50
        err = abs(mp.mpf(P)/Q - mp.mpf(p)/q)
        out.append(float(-(1 + mp.log(err)/mp.log(q))) if err != 0 else float('inf')); mp.mp.dps = old
    return out
