from fractions import Fraction as F
import math
from sympy import prevprime, nextprime
exec(open('consts.py').read().split('def J(u)')[0])   # al, lam, H, Gam(x)
def gamma_in(K,p,M):
    N=3*K//40; h=37*K//40; L0=4*M+10; m=(p-1)//2
    def cnt(a,A): return (A-a)//p+1 if a<=A else 0
    lK=[cnt(a,K)+cnt(p-a,K) for a in range(m+1)]; lN=[cnt(a,N)+cnt(p-a,N) for a in range(m+1)]
    mN=N//p; mK=K//p
    Q=h-L0+3*(N-mN); T=Q//m; E=Q-m*T
    order=sorted(range(1,m+1),key=lambda a:-lK[a])
    eps=[0]*(m+1)
    for a in order[:E]: eps[a]=1
    tot=F(0); Zmin=None
    for a in range(1,m+1):
        b=3*lN[a]; L=T-b+eps[a]; assert L>=0
        tot+=L*(b-F(lK[a]+4,2))+F(L*(L-1),2)
        z=T+eps[a]-F(lK[a]+4,2); Zmin=z if Zmin is None else min(Zmin,z)
    for i in range(L0): tot+=min(2*i+6*mN-mK+F(1,2),Zmin)
    return 2*tot
def Gm(x): return Gam(F(x))
M=40; K=40*8000
for p in (8009, 9001, 12007, 20011, 30011, 50021, 80021, 106661):
    if not (F(K,M)<p<=F(K,3)): continue
    g=gamma_in(K,p,M); pg=p*Gm(F(K,p))
    print(f'K={K} p={p} x=K/p={K/p:.3f}: gamma_in={float(g):.1f}  p*Gamma(K/p)={float(pg):.1f}  diff={float(g-pg):.1f}')
