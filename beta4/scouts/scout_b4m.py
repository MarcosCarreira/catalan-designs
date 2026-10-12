"""Route 3: arcsine array with log kernel, L(i,j) = int dmu t^i int_0^1 s^j (-log s)/(1-ts) ds = sum_u m_{i+u}/(u+j+1)^2,
mu = |t| dt/sqrt(1-t^2). PSLQ of L(i,j) for even i, j against weight <= 3 constants. usage: python3 scout_b4m.py dps"""
import mpmath as m, sys, time
from fractions import Fraction as F
T0=time.time(); m.mp.dps=int(sys.argv[1])
pi,L2,G,z3=m.pi,m.log(2),+m.catalan,m.zeta(3)
I3=m.im(m.polylog(3,m.mpc(1,1)/2)); X3=8*I3-pi**3/16-pi/4*L2**2
B={'1':1,'G':G,'X3':X3,'pi^2':pi**2,'zeta3':z3,'G*log2':G*L2,'pi*log2':pi*L2,'pi^2*log2':pi**2*L2,'log2':L2,'pi':pi,'pi*G':pi*G,'log2^2':L2**2,'log2^3':L2**3}
def h(j,t):
    if abs(t)<m.mpf('0.3'): return m.nsum(lambda u: t**u/(u+j+1)**2,[0,m.inf])
    return (m.polylog(2,t)-sum(t**k/m.mpf(k)**2 for k in range(1,j+1)))/t**(j+1)
def L(i,j):
    return m.re(m.quad(lambda th: m.sin(th)*(m.sin(th)**i*h(j,m.sin(th))+(-m.sin(th))**i*h(j,-m.sin(th))),[0,pi/4,pi/2]))
for (i,j) in [(0,0),(2,0),(4,0),(0,2),(0,4),(2,2),(4,2),(2,4),(6,0),(0,6)]:
    v=L(i,j); ns=list(B); rel=m.pslq([v]+[B[n] for n in ns],maxcoeff=10**8,maxsteps=10**6)
    print((i,j), m.nstr(v,20), None if rel is None else {n:str(F(-c,rel[0])) for c,n in zip(rel[1:],ns) if c}, f"[{time.time()-T0:.0f}s]",flush=True)
