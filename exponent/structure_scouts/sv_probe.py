import sys, time, mpmath as mpm, full_cat
n,b,q,g,h=map(int,sys.argv[1].split(','))
dps=int(sys.argv[2]) if len(sys.argv)>2 else 1800
t=time.time()
Frat,FG=full_cat.build(n,b,q,g,h)
mpm.mp.dps=dps
G=+mpm.catalan
A=mpm.matrix(n,n); B=mpm.matrix(n,n)
for r in range(n):
    for k in range(n):
        x,y=Frat[r][k],FG[r][k]
        A[r,k]=mpm.mpf(x.numerator)/x.denominator+G*mpm.mpf(y.numerator)/y.denominator
        B[r,k]=mpm.mpf(y.numerator)/y.denominator
print('built',round(time.time()-t,1),flush=True)
S=mpm.svd_r(A,compute_uv=False)
sv=sorted([float(mpm.log(s)) for s in S],reverse=True)
print('log sing values of A (desc):',[round(v,1) for v in sv])
d=[sv[i]-sv[i+1] for i in range(len(sv)-1)]
print('gaps:',[round(v,2) for v in d])
print('logdet',round(sum(sv),1),'time',round(time.time()-t,1))
