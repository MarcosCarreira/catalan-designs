import sys, math, time
import numpy as np
from sympy import primerange
def bound44(n,p):
    r=(3*n)%p; first=-(9*n*n-r*r-p*(3*n-r))//p
    INF=10**12; dp=np.full(3*n+1,INF,dtype=np.int64); dp[0]=0
    f=lambda x:x//p
    for rho in range(p):
        M=2+2*f(6*n-rho); c=2*f(6*n)-2*f(6*n-rho)+3*f(2*n+rho)-3*f(2*n)-5
        ind=1 if p>max(6*n-rho,2*n+rho) else 0
        nd=np.full(3*n+1,INF,dtype=np.int64)
        for m in range(0,min(M,3*n)+1):
            Fm=m*m+c*m+3*max(0,m-2*f(6*n-rho))+(ind if m==1 else 0)
            cand=dp[:3*n+1-m]+Fm
            np.minimum(nd[m:],cand,out=nd[m:])
        dp=nd
    return first+int(dp[3*n])
for n in map(int,sys.argv[1].split(',')):
    T=time.time(); tot=0; small=0
    for p in primerange(2,6*n+1):
        if p*p>8*n: b=max(bound44(n,p),-12*n)
        else: b=-12*n*int(math.floor(math.log(8*n)/math.log(p)+1e-12)); small+=b*math.log(p)
        tot+=b*math.log(p)
    print(f'n={n}: sum of Lemma 4.1/4.4 bounds * log p / n^2 = {tot/n**2:.3f} (of which p <= sqrt(8n): {small/n**2:.3f}) [{time.time()-T:.1f}s]',flush=True)
