"""Route 3: even weight |sin th| (th^2 - pi |th|) (moments rational: the pi-part of the th^2 moments equals p_n).
Moments check, then starting values S0, S1, P1 of the kernel 1/(1-ts) as in scout_b4f (family A with this weight).
usage: python3 scout_b4j.py <mom|S0|S1|P1> <dps>"""
import mpmath as m, sys, time
from fractions import Fraction as F
T0=time.time(); which=sys.argv[1]; m.mp.dps=int(sys.argv[2])
pi,L2,G,z3=m.pi,m.log(2),m.catalan,m.zeta(3)
q=lambda t: t*t-pi*t       # on (0, pi/2); weight is even in theta
if which=='mom':
    for n in range(1,12,2):
        v=m.quad(lambda t: q(t)*m.sin(t)**n,[0,pi/2])
        print(n, m.nstr(v,20), m.pslq([v,1,pi,pi**2],maxcoeff=10**12,maxsteps=10**6), m.identify(v))
    sys.exit()
b4=(m.zeta(4,m.mpf(1)/4)-m.zeta(4,m.mpf(3)/4))/256
wq=m.mpc(1,1)/2; I4=m.im(m.polylog(4,wq)); I3=m.im(m.polylog(3,wq))
B={'beta4':b4,'pi^2*G':pi**2*G,'G^2':G**2,'G*log2^2':G*L2**2,'pi*G*log2':pi*G*L2,'ImLi4':I4,'log2*ImLi3':L2*I3,'pi*ImLi3':pi*I3,
   'pi^4':pi**4,'pi^2*log2^2':pi**2*L2**2,'log2^4':L2**4,'pi^3*log2':pi**3*L2,'pi*log2^3':pi*L2**3,'pi*zeta3':pi*z3,'zeta3*log2':z3*L2,
   'Li4half':m.polylog(4,m.mpf(1)/2)}
B3={'pi^3':pi**3,'pi*log2^2':pi*L2**2,'G*log2':G*L2,'ImLi3':I3,'pi*G':pi*G,'zeta3':z3,'pi^2*log2':pi**2*L2,'log2^3':L2**3}
B1={'1':1,'pi':pi,'log2':L2,'pi^2':pi**2,'pi*log2':pi*L2,'log2^2':L2**2,'G':G}
if which in ('P1','S1'): B={**B,**B3,**B1}
at=lambda t: m.log(m.tan(t/2+pi/4))
s=m.sin
f={'S0':lambda t: s(t)*q(t)*2*at(t)/s(t),                       # Li1(y)/y + Li1(-y)/(-y) = 2 artanh(y)/y
   'S1':lambda t: s(t)*q(t)*(-2*m.log(m.cos(t))),                # Li1(y) + Li1(-y) = -log(1-y^2)
   'P1':lambda t: s(t)*q(t)*(-2*m.log(m.cos(t)))/s(t)**2}[which]  # (Li1(y)-y)/y^2 + (Li1(-y)+y)/y^2
v=m.re(m.quad(f,[0,pi/8,pi/4,3*pi/8,pi/2]))
ns=list(B); rel=m.pslq([v]+[B[n] for n in ns],maxcoeff=10**6,maxsteps=10**6)
print(f"{which} dps={m.mp.dps}: {m.nstr(v,25)} ->", None if rel is None else {n:str(F(-c,rel[0])) for c,n in zip(rel[1:],ns) if c}, f"[{time.time()-T0:.0f}s]")
