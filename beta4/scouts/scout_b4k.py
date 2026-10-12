"""Route 3: boundary sequences of the measure nu = |t| (arcsin^2 t - pi |arcsin t|) dt/sqrt(1-t^2):
K-_i = int t^(i-1) (-log(1-t)) dnu,  K+_j = int int_0^1 s^j/(1-ts) ds dnu  (K-_0 = K+_0 = S0, K-_1 = S1, K+_1 = P1).
PSLQ against {1, beta4, G, G^2, pi^2, pi^4}. usage: python3 scout_b4k.py <dps> <imax>"""
import mpmath as m, sys, time, json
from fractions import Fraction as F
T0=time.time(); m.mp.dps=int(sys.argv[1]); IMAX=int(sys.argv[2])
pi,G=m.pi,m.catalan
b4=(m.zeta(4,m.mpf(1)/4)-m.zeta(4,m.mpf(3)/4))/256
B={'1':1,'beta4':b4,'G':G,'G^2':G**2,'pi^2':pi**2,'pi^4':pi**4}
q=lambda th: th*th-pi*th
def lm(th): return -m.log(2)-2*m.log(m.sin(pi/4-th/2))   # -log(1-sin th), stable at pi/2
def lp(th): return -m.log(1+m.sin(th))
def nu_int(Fp):   # Fp(th) = F(s)+F(-s), s = sin th;  int F dnu = int_0^{pi/2} sin(th) q(th) Fp(th) dth
    return m.re(m.quad(lambda th: m.sin(th)*q(th)*Fp(th),[0,pi/8,pi/4,3*pi/8,pi/2]))
def Km(i):
    if i==0: return lambda th: (lm(th)-lp(th))/m.sin(th)
    return lambda th: m.sin(th)**(i-1)*lm(th)+(-m.sin(th))**(i-1)*lp(th)
def gpair(j,th):
    s_=m.sin(th); tot=0
    for t,L in ((s_,lm(th)),(-s_,lp(th))):
        if abs(t)<m.mpf('0.25'): tot+=m.nsum(lambda u: t**u/(u+j+1),[0,m.inf])
        else: tot+=(L-sum(t**k/k for k in range(1,j+1)))/t**(j+1)
    return tot
def Kp(j): return lambda th: gpair(j,th)
G=+G
B={'1':1,'beta4':b4,'G':G,'G^2':G**2,'pi^2':pi**2,'pi^4':pi**4}
out={}
for nm,fn in [('K-',Km),('K+',Kp)]:
    for i in range(IMAX+1):
        v=nu_int(fn(i)); ns=list(B)
        if abs(v)<m.mpf(10)**(-m.mp.dps+15): print(f'{nm}_{i} = 0',flush=True); out[f'{nm}{i}']='0'; continue
        rel=m.pslq([v]+[B[n] for n in ns],maxcoeff=10**9,maxsteps=10**6)
        d=None if rel is None else {n:str(F(-c,rel[0])) for c,n in zip(rel[1:],ns) if c}
        out[f"{nm}{i}"]=d
        print(f"{nm}_{i} = {m.nstr(v,20)}: {d}  [{time.time()-T0:.0f}s]",flush=True)
json.dump(out,open(f'b4k_{m.mp.dps}.json','w'),indent=1)
