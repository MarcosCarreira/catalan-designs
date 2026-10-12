import sys, math
sys.argv=['x']
src=open('fz5_check.py').read().split('if __name__')[0]
g={'__file__':'fz5_check.py','__name__':'m'}; exec(src,g)
from flint import fmpq, fmpq_poly, fmpz_mat
from sympy import primerange
def vp(x,p):
    x=int(x); v=0
    if x==0: return 10**9
    while x%p==0: x//=p; v+=1
    return v
def vq(q,p): return vp(q.p,p)-vp(q.q,p)
n=int(sys.argv[1]) if len(sys.argv)>1 else 1
K,N,h,a,b,lead=g['entries'](n)
L=1
for q in a+b: L=L*int(q.q)//math.gcd(L,int(q.q))
ai=[int((q*L).p) for q in a]; bi=[int((q*L).p) for q in b]
xs=list(range(h+1)); ys=[fmpq(int(fmpz_mat(h,h,[ai[i+j]+x*bi[i+j] for i in range(h) for j in range(h)]).det())) for x in xs]
coef=ys[:]
for k in range(1,h+1):
    for i in range(h,k-1,-1): coef[i]=(coef[i]-coef[i-1])/(xs[i]-xs[i-k])
t=fmpq_poly([0,1]); P=fmpq_poly([coef[h]])
for k in range(h-1,-1,-1): P=P*(t-xs[k])+coef[k]
Delta=P/fmpq(L)**h
fac=math.factorial
SK=fmpq(fac(K)**(2*h)*4**(h-1), fac(N)**(12*h)*math.prod(fac(2*i)**2 for i in range(1,h)))
F=Delta*SK
def vG(Pol,p): return min(vq(Pol[k],p) for k in range(Pol.degree()+1) if Pol[k]!=0)
print(f"n={n} K={K} N={N} h={h}")
bad=0
for p in primerange(2,3*K):
    vF=vG(F,p); bound=-6*h*int(math.floor(math.log(5*K)/math.log(p)+1e-12))-h*vp(24,p)
    ok312= vF>=bound
    line=f"p={p:3d} vG(F_K)={vF:6d}  (3.12) bound {bound:6d} {'ok' if ok312 else 'FAIL'}"
    if not ok312: bad+=1
    vD=vG(Delta,p)
    if p>K:
        line+=f" | p>K: vG(Delta)={vD} {'ok' if vD>=0 else 'FAIL'}"; bad+= vD<0
    elif p>=7 and p<=K<3*p and p*p>2*K and 2*N<p and 5*N<=2*p-2:
        v=K-p*(K//p); u=max(0,N+v-p+1); tp=min(N,v)+u; rp=max(0,K+4*N-2*p+2)
        g_out=(-7*(K-p)+6*tp-1-min(rp,p-1-N+u)) if K<2*p else (-7*(K-p)+3+12*N+5*tp-min(rp,p+u))
        line+=f" | (4.9) holds: vG(Delta)={vD} gamma_out={g_out} {'ok' if vD>=g_out else 'FAIL'}"; bad+= vD<g_out
    print(line)
print('violations:',bad)
