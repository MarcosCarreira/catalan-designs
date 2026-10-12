import mpmath as mp
from nu_arrays import *
mp.mp.dps = 50
pi = mp.pi; G = +mp.catalan; b4 = (mp.zeta(4, mp.mpf(1)/4)-mp.zeta(4, mp.mpf(3)/4))/256
val = {'1': 1, 'G': G, 'b4': b4}
m, km, nu, K = build_nu(14)
q = lambda th: th*th - pi*th
# moments check
for k in range(4):
    v = 2*mp.quad(lambda th: q(th)*mp.sin(th)**(2*k+1), [0, pi/2]); print('nu', 2*k, mp.nstr(v, 20), nu[2*k])
def lm(th): return -mp.log(2)-2*mp.log(mp.sin(pi/4-th/2))
def lp(th): return -mp.log(1+mp.sin(th))
def gfun(j, t, L):
    if abs(t) < mp.mpf('0.25'): return mp.nsum(lambda u: t**u/(u+j+1), [0, mp.inf])
    return (L - sum(t**k/k for k in range(1, j+1)))/t**(j+1)
def Mnum(i, j):
    f = lambda th: mp.sin(th)*q(th)*(mp.sin(th)**i*gfun(j, mp.sin(th), lm(th)) + (-mp.sin(th))**i*gfun(j, -mp.sin(th), lp(th)))
    return mp.re(mp.quad(f, [0, pi/8, pi/4, 3*pi/8, pi/2]))
worst = 0
for (i, j) in [(0,0),(2,0),(3,0),(5,0),(7,0),(8,0),(4,2),(6,2),(9,4),(10,6),(12,4),(13,8),(14,14)]:
    ex = M_nu(i, j, nu, K); e = sum(mp.mpf(v.numerator)/v.denominator*val[k] for k, v in ex.items())
    n = Mnum(i, j); worst = max(worst, abs(e-n)); print((i, j), mp.nstr(n, 25), mp.nstr(e-n, 3), {k: str(v) for k, v in ex.items() if k != '1'})
print('worst', mp.nstr(worst, 3))
