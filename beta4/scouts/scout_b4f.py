import mpmath as m, time, json
from fractions import Fraction as F
T0=time.time()
m.mp.dps=130
pi,L2,G,z3=m.pi,m.log(2),m.catalan,m.zeta(3)
b4=(m.zeta(4,m.mpf(1)/4)-m.zeta(4,m.mpf(3)/4))/256
I4=m.im(m.polylog(4,m.mpc(0.5,0.5))); I3=m.im(m.polylog(3,m.mpc(0.5,0.5))); Li4h=m.polylog(4,m.mpf(1)/2)
EPS=m.mpf(10)**-125
def Li(k,y):
    if k==1: return -m.log(1-y) if y<0.5 else -(m.log(2)+2*m.log(m.cos(m.asin(y)/2+pi/4)))  # stable near y=1
    return m.re(m.polylog(k,y))
def sym(f,wt):   # integral over (-pi/2, pi/2) of wt(theta) f(sin theta)
    return m.re(m.quad(lambda th: m.re(wt(th)*f(m.sin(th)) + wt(-th)*f(-m.sin(th))),[0,pi/4,pi/2-EPS]))
W={'A':lambda th: abs(m.sin(th)), 'B':lambda th: th}
def starts(fam,r):
    w=W[fam]
    S0=sym(lambda y: Li(r+1,y)/y, w)
    S1=sym(lambda y: Li(r+1,y), w)
    P1=sym(lambda y: (Li(r+1,y)-y)/(y*y), w)
    return {'S0':S0,'S1':S1,'P1':P1}
basis={'1':1,'pi':pi,'log2':L2,'G':G,'pi^2':pi**2,'pi*log2':pi*L2,'log2^2':L2**2,
 'pi^3':pi**3,'pi*log2^2':pi*L2**2,'G*log2':G*L2,'ImLi3':I3,'zeta3':z3,'pi^2*log2':pi**2*L2,'log2^3':L2**3,'pi*G':pi*G,
 'beta4':b4,'ImLi4':I4,'ImLi3*log2':I3*L2,'G*log2^2':G*L2**2,'pi^2*G':pi**2*G,'pi^3*log2':pi**3*L2,'pi*log2^3':pi*L2**3,'pi*zeta3':pi*z3,
 'pi^4':pi**4,'zeta3*log2':z3*L2,'pi^2*log2^2':pi**2*L2**2,'log2^4':L2**4,'Li4half':Li4h}
names=list(basis)
out={}
import sys
FAMS={'A':(0,1,2),'B':(0,1)}
sel=sys.argv[1]
for fam,rs in [(sel[0],tuple(int(c) for c in sel[1:]))]:
    for r in rs:
        st=starts(fam,r)
        for k,v in st.items():
            rel=m.pslq([v]+[basis[n] for n in names],maxcoeff=10**9,maxsteps=10**7)
            if rel is None: print(f"{fam} r={r} {k}: NO RELATION",flush=True); continue
            vec={n:str(F(-c,rel[0])) for c,n in zip(rel[1:],names) if c}
            out[f"{fam}{r}{k}"]=vec
            print(f"{fam} r={r} {k} = {m.nstr(v,12)}: {vec}  [{time.time()-T0:.0f}s]",flush=True)
json.dump(out,open(f'starts_{sel}.json','w'),indent=1)
