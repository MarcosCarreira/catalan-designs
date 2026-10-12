import sys, math, time
from flint import fmpq, fmpq_mat
import mpmath as mp
from functools import lru_cache
@lru_cache(None)
def H(k,s): return fmpq(0) if k==0 else H(k-1,s)+fmpq(1,k**s)
def Lv(n,a,b,c):
    N6=6*n; rows=[]
    for k in range(N6+1):
        w=fmpq(math.comb(N6,k)**2*math.comb(2*n+k,k)**3); A=3*a*H(k,2)+2*b*H(k,3)
        rows.append((w,-A,(3*H(2*n+k,1)+2*H(N6-k,1)-5*H(k,1))*A,-a*H(k,3)-b*H(k,4)))
    out=[]
    for s in range(6*n-1):
        tot=fmpq(0) if not isinstance(c,mp.mpf) else mp.mpf(0)
        for k,(w,U,V,h) in enumerate(rows):
            x=-k; P=x**s; dP=s*x**(s-1) if s>0 else 0
            if isinstance(c,mp.mpf):
                tot+=mp.mpf(int((w*(U*dP+V*P+6*h*P)).p))/int((w*(U*dP+V*P+6*h*P)).q)+6*int(w.p)*P*c
            else: tot+=w*(U*dP+V*P+6*(c+h)*P)
        out.append(tot)
    return out
def logden(n,a,b,c):
    L=Lv(n,a,b,fmpq(c)); m=3*n
    d=fmpq_mat(m,m,[L[i+j] for i in range(m) for j in range(m)]).det()
    norm=math.prod(math.factorial(i)**2 for i in range(m))
    J=d/norm; return math.log(int(J.q)), J
mode=sys.argv[1]
if mode=='den':
    for n in map(int,sys.argv[2].split(',')):
        T=time.time(); ld,J=logden(n,1,1,1); print(f'n={n}: -log den(J_n)/n^2 = {-ld/n**2:.3f}  [{time.time()-T:.1f}s]',flush=True)
else:
    for n in map(int,sys.argv[2].split(',')):
        T=time.time(); mp.mp.dps=200+40*n*n
        c=mp.zeta(3)+mp.zeta(4); L=Lv(n,1,1,c); m=3*n
        M=mp.matrix(m,m)
        for i in range(m):
            for j in range(m): M[i,j]=L[i+j]
        J=mp.det(M)/math.prod(math.factorial(i)**2 for i in range(m))
        print(f'n={n}: real c = zeta3+zeta4: log|J_n|/n^2 = {mp.nstr(mp.log(abs(J))/n**2,6)}  [{time.time()-T:.1f}s]',flush=True)
