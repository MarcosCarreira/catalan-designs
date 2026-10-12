"""Route 3: starting values S0, S1, P1 (kernel 1/(1-ts), as in scout_b4f) for odd weights w(theta), with
homogeneous bases; usage: python3 scout_b4i.py <weight: 1|3> <S0|S1|P1> <dps>"""
import mpmath as m, sys, time
T0=time.time()
k=int(sys.argv[1]); which=sys.argv[2]; m.mp.dps=int(sys.argv[3])
pi,L2,G,z3=m.pi,m.log(2),m.catalan,m.zeta(3)
b4=(m.zeta(4,m.mpf(1)/4)-m.zeta(4,m.mpf(3)/4))/256
wq=m.mpc(1,1)/2; I4=m.im(m.polylog(4,wq)); I3=m.im(m.polylog(3,wq))
if k==1:
    W=lambda t: t
    B={'pi^3':pi**3,'pi*log2^2':pi*L2**2,'G*log2':G*L2,'ImLi3':I3,'pi*G':pi*G,'zeta3':z3,'pi^2*log2':pi**2*L2,'log2^3':L2**3}
else:
    W=lambda t: t**3-3*pi**2/4*t
    B={'beta4':b4,'pi^2*G':pi**2*G,'G^2':G**2,'G*log2^2':G*L2**2,'pi*G*log2':pi*G*L2,'ImLi4':I4,'log2*ImLi3':L2*I3,'pi*ImLi3':pi*I3,
       'pi^4':pi**4,'pi^2*log2^2':pi**2*L2**2,'log2^4':L2**4,'pi^3*log2':pi**3*L2,'pi*log2^3':pi*L2**3,'pi*zeta3':pi*z3,'zeta3*log2':z3*L2,
       'Li4half':m.polylog(4,m.mpf(1)/2)}
    if which=='P1':
        B={'pi^3':pi**3,'pi*log2^2':pi*L2**2,'G*log2':G*L2,'ImLi3':I3,'pi*G':pi*G,'zeta3':z3,'pi^2*log2':pi**2*L2,'log2^3':L2**3}
at=lambda t: m.log(m.tan(t/2+pi/4))
f={'S0':lambda t: -2*W(t)*m.log(m.cos(t))/m.sin(t),
   'S1':lambda t: 2*W(t)*at(t),
   'P1':lambda t: 2*W(t)*(at(t)-m.sin(t))/m.sin(t)**2}[which]
v=m.re(m.quad(f,[0,pi/8,pi/4,3*pi/8,pi/2]))
ns=list(B); rel=m.pslq([v]+[B[n] for n in ns],maxcoeff=10**6,maxsteps=10**6)
from fractions import Fraction as F
print(f"w=theta^{k} {which} dps={m.mp.dps}: {m.nstr(v,25)} ->",
      None if rel is None else {n:str(F(-c,rel[0])) for c,n in zip(rel[1:],ns) if c}, f"[{time.time()-T0:.0f}s]")
