"""H1 calibration: Rivoal-Zudilin q = 6 forms S_n(-1) = c2 beta(2) + c4 beta(4) + c0 (exact), elimination of beta(2)
between two parameter families, denominators versus size."""
import sys, time, math
from fractions import Fraction as F
import mpmath as m
exec(open('rz_scout.py').read().split("etas=")[0])   # pmul, shift, series_div, forms (P_j)
T0=time.time()
m.mp.dps=60
b2=m.catalan; b4=(m.zeta(4,m.mpf(1)/4)-m.zeta(4,m.mpf(3)/4))/256
def full_form(eta,n):
    """returns (c0, c2, c4) with S_n(-1) = c2*beta(2) + c4*beta(4) + c0, using eqs. (12), (15), (16)."""
    q=len(eta)-1; e0=eta[0]; h0=e0*n+1; hq=eta[q]*n+1
    # recompute coefficients c_kj (copy of forms() internals, returning coef)
    h=[None]+[F(2*eta[j]*n+1,2) for j in range(1,q)]+[hq]
    gam=F(1)
    for j in range(1,q): gam*=math.factorial((e0-2*eta[j])*n)
    gam/=math.factorial(eta[q]*n)**2
    num=[F(h0)*gam, 2*gam]
    for i in range(hq-1): num=pmul(num,[F(1+i),F(1)])
    for i in range(hq-1): num=pmul(num,[F(1+h0-hq+i),F(1)])
    poles={}
    for j in range(1,q):
        L=1+h0-int(2*h[j])
        for i in range(L):
            p=-(h[j]+i); poles[p]=poles.get(p,0)+1
    roots=[]
    for p,ml in poles.items(): roots+=[p]*ml
    mm=q-1; coef={}
    for p,ml in poles.items():
        den=[F(1)]
        for r in roots:
            if r!=p: den=pmul(den,[-r,F(1)])
        nu=[F(0)]*(mm-ml)+shift(num,p); de=shift(den,p)
        ser=series_div(nu,de,mm); k=-p+F(1,2)
        for j in range(1,mm+1): coef[(k,j)]=ser[mm-j]
    # sum_{t>=0} (-1)^t/(t+k-1/2)^j = (-1)^(k+1) 2^j beta(j) - (-1)^k sum_{l=1}^{k-1} (-1)^l/(l-1/2)^j
    cb={j: F(0) for j in range(1,q)}; c0=F(0)
    for (k,j),c in coef.items():
        k=int(k); sg=(-1)**(k+1)
        cb[j]+=c*sg*2**j
        c0-=c*(-1)**k*sum(F((-1)**l)/F(2*l-1,2)**j for l in range(1,k))
    assert cb[1]==0 and cb[3]==0 and cb[5]==0, cb
    return (c0, cb[2], cb[4]), coef
def direct(eta,n):
    """S_n(-1) by summing the series (6) numerically (positive control)."""
    q=len(eta)-1; e0=eta[0]; h0=e0*n+1
    h=[h0]+[m.mpf(2*eta[j]*n+1)/2 for j in range(1,q)]+[eta[q]*n+1]
    g=m.mpf(1)
    for j in range(1,q): g*=m.factorial((e0-2*eta[j])*n)
    g/=m.factorial(eta[q]*n)**2
    def term(t):
        v=(2*t+h0)*m.gamma(t+h0)/m.gamma(t+1)
        for j in range(1,q+1): v*=m.gamma(t+h[j])/m.gamma(t+1+h0-h[j])
        return g*v*(-1)**int(t)
    return m.nsum(term,[0,m.inf])
def den(x): return x.denominator
ETA1=(3,1,1,1,1,1,1); ETA2=(5,2,2,2,2,2,2)
print("positive control: exact form vs direct series")
for eta in (ETA1,ETA2):
    for n in (1,2):
        (c0,c2,c4),_=full_form(eta,n)
        ex=c2*b2+c4*b4+c0
        print(f"  eta={eta} n={n}: form {m.nstr(m.mpf(c2.numerator)/c2.denominator*b2+m.mpf(c4.numerator)/c4.denominator*b4+m.mpf(c0.numerator)/c0.denominator,20)}  direct {m.nstr(direct(eta,n),20)}",flush=True)
print("H1: L_n = c2(eta2) S_n(eta1) - c2(eta1) S_n(eta2) = A_n beta(4) + B_n")
for n in range(1,int(sys.argv[1]) if len(sys.argv)>1 else 7):
    (a0,a2,a4),_=full_form(ETA1,n); (b0,bb2,bb4),_=full_form(ETA2,n)
    A=bb2*a4-a2*bb4; B=bb2*a0-a2*b0
    D=A.denominator*B.denominator//math.gcd(A.denominator,B.denominator)
    Lval=m.mpf(A.numerator)/A.denominator*b4+m.mpf(B.numerator)/B.denominator
    sz=lambda x: m.log(abs(m.mpf(x.numerator)/x.denominator)) if x!=0 else -m.inf
    print(f"  n={n}: log|L_n|/n = {m.nstr(m.log(abs(Lval))/n,6)}, log D_n/n = {m.nstr(m.log(D)/n,6)}, sum = {m.nstr((m.log(abs(Lval))+m.log(D))/n,6)}; "
          f"log|S1|/n = {m.nstr(m.log(abs(a2*b2+a4*b4+a0 if False else m.mpf(1)))/n,3)}  [{time.time()-T0:.0f}s]",flush=True)
