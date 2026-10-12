from fractions import Fraction as F
import math
al=F(3,40); lam=F(37,40); H=1+2*al
fl=math.floor
def pos(v): return v if v>0 else 0
def zint(f, x):
    """integral over z in (0,1/2) of f(z), f piecewise constant with breaks where frac(x)+-z or frac(al x)+-z cross integers"""
    br={F(0),F(1,2)}
    for u in (x, al*x):
        fr=u-fl(u)
        for c in (fr, 1-fr):
            if 0<c<F(1,2): br.add(c)
    br=sorted(br); tot=F(0)
    for a,b in zip(br,br[1:]):
        if b>a: tot+=(b-a)*f((a+b)/2)
    return tot
def ell(u,z): return fl(u-z)+fl(u+z)+1
def Gam(x):
    T=fl(2*H*x); s=H*x-F(T,2); q=fl(2*x); npl=(2*x-q)/2
    I=zint(lambda z:(T-3*ell(al*x,z))*(T+3*ell(al*x,z)-ell(x,z)-5), x)
    return I+s*(2*T-q-5)+pos(s-npl)
def J(u):
    m=fl(2*u); return m*u-F(m*(m+1),4)
def Nf(x): return 2*lam*x*fl(x)-12*lam*x*fl(al*x)-2*J(lam*x)
def R(x): return -Gam(x)-Nf(x)
# integral of R(x)/x^3 over [3,20]: R is piecewise polynomial of degree <= 2 in x between breakpoints; use 3-point Gauss on fine pieces
import mpmath as mp
mp.mp.dps=30
br=set()
for den in (F(1,1),F(1,2),1/(2*H),1/(2*lam),1/al,1/(2*al)):
    k=0
    while k*den<=20:
        if k*den>=3: br.add(k*den)
        k+=1
br.add(F(3)); br.add(F(20)); br=sorted(br)
# add kink points where s = n+ (s - n+ changes sign): solve inside each piece by bisection-free refinement: subdivide each piece into 64
gx,gw=mp.gauss_legendre if False else (None,None)
nodes=[(-mp.sqrt(mp.mpf(3)/5),mp.mpf(5)/9),(0,mp.mpf(8)/9),(mp.sqrt(mp.mpf(3)/5),mp.mpf(5)/9)]
def integ(a,b,sub=64):
    tot=mp.mpf(0)
    for i in range(sub):
        lo=a+(b-a)*i/sub; hi=a+(b-a)*(i+1)/sub; m=(lo+hi)/2; h=(hi-lo)/2
        for t,w in nodes:
            x=F(m+h*t) if isinstance(m,F) else None
            xx=m+h*F(str(t)) if False else None
            xv=mp.mpf(m.numerator)/m.denominator+mp.mpf(h.numerator)/h.denominator*t
            xq=F(str(mp.nstr(xv,28)))
            tot+=w*mp.mpf(h.numerator)/h.denominator*(mp.mpf(R(xq).numerator)/R(xq).denominator)/xv**3
    return tot
I320=sum(integ(a,b) for a,b in zip(br,br[1:]))
print('int_3^20 R/x^3 =',mp.nstr(I320,12),' paper:',float(F(322437603634266857629,7535670527041937280000)))
