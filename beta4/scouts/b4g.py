import mpmath as m, time, sys
T0=time.time()
m.mp.dps=int(sys.argv[2]) if len(sys.argv)>2 else 160
pi,L2,G,z3=m.pi,m.log(2),m.catalan,m.zeta(3)
b4=(m.zeta(4,m.mpf(1)/4)-m.zeta(4,m.mpf(3)/4))/256
I4=m.im(m.polylog(4,m.mpc(0.5,0.5))); I3=m.im(m.polylog(3,m.mpc(0.5,0.5)))
def Li(k,th,sgn):   # Li_k(sgn*sin th), stable near sin th -> 1
    y=sgn*m.sin(th)
    if k==1:
        if sgn>0: return -(L2+2*m.log(m.cos(th/2+pi/4)))
        return -m.log(1-y)
    return m.re(m.polylog(k,y))
def B_int(f):  # int_{-pi/2}^{pi/2} theta * f(theta), f given for +theta and -theta
    return m.quad(lambda th: th*(f(th,1)-f(th,-1)),[0,pi/8,pi/4,3*pi/8,pi/2])
which=sys.argv[1]
r=int(which[1]); kind=which[2:]
if kind=='S0': v=B_int(lambda th,s: Li(r+1,th,s)/(s*m.sin(th)))
if kind=='S1': v=B_int(lambda th,s: Li(r+1,th,s))
if kind=='P1': v=B_int(lambda th,s: (Li(r+1,th,s)-s*m.sin(th))/(m.sin(th)**2))
odd={'pi':pi,'G':G,'pi*log2':pi*L2,'pi^3':pi**3,'pi*log2^2':pi*L2**2,'G*log2':G*L2,'ImLi3':I3,
     'beta4':b4,'ImLi4':I4,'ImLi3*log2':I3*L2,'G*log2^2':G*L2**2,'pi^2*G':pi**2*G,'pi^3*log2':pi**3*L2,'pi*log2^3':pi*L2**3,'pi*zeta3':pi*z3}
ev={'1':1,'log2':L2,'pi^2':pi**2,'log2^2':L2**2,'zeta3':z3,'pi^2*log2':pi**2*L2,'log2^3':L2**3,'pi^4':pi**4,'zeta3*log2':z3*L2,'pi^2*log2^2':pi**2*L2**2,'log2^4':L2**4,'Li4half':m.polylog(4,0.5)}
from fractions import Fraction as F
for nm,B in (('odd',odd),('even',ev)):
    ns=list(B); rel=m.pslq([v]+[B[n] for n in ns],maxcoeff=10**6,maxsteps=10**7)
    print(f"B {which} = {m.nstr(v,15)} [{nm} basis]:", None if rel is None else {n:str(F(-c,rel[0])) for c,n in zip(rel[1:],ns) if c}, f"[{time.time()-T0:.0f}s]",flush=True)
