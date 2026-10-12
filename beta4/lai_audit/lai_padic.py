import sys, math, time
from sympy import primerange
exec(open('lai_trend.py').read().split("mode=sys.argv[1]")[0])
def vp(x,p):
    x=abs(int(x)); v=0
    if x==0: return 10**9
    while x%p==0: x//=p; v+=1
    return v
def F(rho,m,n,p):
    f=lambda x: x//p
    c=2*f(6*n)-2*f(6*n-rho)+3*f(2*n+rho)-3*f(2*n)-5
    return m*m+c*m+3*max(0,m-2*f(6*n-rho))+(1 if (m==1 and p>max(6*n-rho,2*n+rho)) else 0)
def bound44(n,p):
    r=(3*n)%p
    first=-(9*n*n-r*r-p*(3*n-r))//p   # exact: -(9n^2/p - 3n + p{}(1-{}))
    INF=10**9; dp=[0]+[INF]*(3*n)
    for rho in range(p):
        M=2+2*((6*n-rho)//p); nd=[INF]*(3*n+1)
        Fv=[F(rho,m,n,p) for m in range(M+1)]
        for s in range(3*n+1):
            if dp[s]>=INF: continue
            for m in range(min(M,3*n-s)+1):
                v=dp[s]+Fv[m]
                if v<nd[s+m]: nd[s+m]=v
        dp=nd
    return first+dp[3*n]
for n in map(int,sys.argv[1].split(',')):
    T=time.time(); ld,J=logden(n,1,1,1)
    tot_meas=0; tot_bound=0; viol=[]; rows=[]
    for p in primerange(2,6*n*2):
        v=vp(J.p,p)-vp(J.q,p)
        if p>6*n: b=0
        elif p*p>8*n: b=bound44(n,p)
        else: b=-12*n*int(math.floor(math.log(8*n)/math.log(p)+1e-12))
        if v<b: viol.append((p,v,b))
        tot_meas+=min(v,0)*math.log(p); tot_bound+=b*math.log(p)
        rows.append((p,v,b))
    tight=sum(1 for p,v,b in rows if v==b and p<=6*n)
    print(f'n={n}: violations {viol}; primes with v_p = bound: {tight}; sum of bounds*log p/n^2 = {tot_bound/n**2:.3f}; measured -log den/n^2 = {-ld/n**2:.3f} [{time.time()-T:.1f}s]',flush=True)
