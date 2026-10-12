"""Rivoal-Zudilin very-well-poised forms S_n(-1) = sum_l P_{2l}(-1) 2^{-2l} beta(2l) - P_0(-1), q = 6 (Rivoal and Zudilin, Math. Ann. 326 (2003); eqs. 4-16 of the version used).
Exact partial fractions of R_n(t) (eq. 7/9); coefficients of beta(2), beta(4) for small n."""
from fractions import Fraction as F
from math import factorial
import sys, itertools
def pmul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b): out[i+j]+=x*y
    return out
def lin(c): return [F(c),F(1)]      # t + c
def shift(a,p):   # coefficients of a(u+p)
    a=list(a); n=len(a)
    for i in range(n):
        for j in range(n-2,i-1,-1): a[j]+=p*a[j+1]
    return a
def series_div(num,den,N):  # power series num/den up to u^{N-1}, den[0]!=0
    out=[]; r=list(num)+[F(0)]*N
    for k in range(N):
        c=r[k]/den[0]; out.append(c)
        for j,d in enumerate(den):
            if k+j<len(r): r[k+j]-=c*d
    return out
def forms(eta,n):
    q=len(eta)-1; e0=eta[0]
    h0=e0*n+1; h=[None]+[F(2*eta[j]*n+1,2) for j in range(1,q)]+[eta[q]*n+1]
    gam=F(1)
    for j in range(1,q): gam*=factorial((e0-2*eta[j])*n)
    gam/=factorial(eta[q]*n)**2
    # numerator: gam (2t+h0) (t+1)_{hq-1} (t+1+h0-hq)_{hq-1}; denominator factors with poles
    hq=h[q]
    num=[gam*2*F(h0,2)*1, gam*2]  # 2t+h0 = 2(t + h0/2)
    num=[F(h0)*gam, 2*gam]
    for i in range(hq-1): num=pmul(num,lin(1+i))
    for i in range(hq-1): num=pmul(num,lin(1+h0-hq+i))
    poles={}   # pole location p (t=p) -> multiplicity
    for j in range(1,q):
        L=1+h0-int(2*h[j])  # length of (t+hj)_{1+h0-2hj}
        for i in range(L):
            p=-(h[j]+i); poles[p]=poles.get(p,0)+1
    den_lin=[]  # list of roots
    for p,mlt in poles.items(): den_lin+=[p]*mlt
    m=q-1
    coef={}  # (k, j) -> c_kj  with pole at t = -k + 1/2
    for p,mlt in poles.items():
        # R*(t-p)^m = num(t) * (t-p)^(m-mlt) / prod_{other roots r}(t-r)
        others=[r for r in den_lin if r!=p]
        den=[F(1)]
        for r in others: den=pmul(den,[-r,F(1)])
        nu=shift(num,p); de=shift(den,p)
        nu=[F(0)]*(m-mlt)+nu
        ser=series_div(nu,de,m)          # Taylor coeffs of g(u) = R*(u)^m at u=0
        k=-p+F(1,2)
        for j in range(1,m+1):
            coef[(k,j)]=ser[m-j]          # coefficient of 1/(t-p)^j
    P={j: sum(c*(-1)**int(k-1) for (k,jj),c in coef.items() if jj==j) for j in range(1,q)}
    return P
etas=[(3,1,1,1,1,1,1),(4,1,1,1,1,1,1),(4,1,2,1,2,1,1),(5,2,2,2,2,2,2),(4,1,1,1,1,1,2)]
for eta in etas:
    if not (sum(eta[1:])<=F(len(eta)-2,2)*eta[0] and all(0<e<F(eta[0],2) for e in eta[1:])): print(eta,"violates (3)"); continue
    rows=[]
    for n in range(1,5):
        P=forms(eta,n)
        rows.append((n, P[1], P[3], P[2]/4, P[4]/16))
    print(eta)
    for r in rows: print(f"   n={r[0]}: P1(-1)={r[1]}, P3(-1)={r[2]}, coef beta(2)={r[3]}, coef beta(4)={r[4]}")
