import sys, time, mpmath as mpm, full_cat, math
for arg in sys.argv[1:]:
    n,b,q,g,h,dps=map(int,arg.split(','))
    t=time.time(); Frat,FG=full_cat.build(n,b,q,g,h); mpm.mp.dps=dps; G=+mpm.catalan
    A=mpm.matrix(n,n); B=mpm.matrix(n,n)
    for r in range(n):
        for k in range(n):
            x,y=Frat[r][k],FG[r][k]; yy=mpm.mpf(y.numerator)/y.denominator
            A[r,k]=mpm.mpf(x.numerator)/x.denominator+G*yy; B[r,k]=yy
    U,S,V=mpm.svd_r(A)          # A = U diag(S) V
    i=min(range(n),key=lambda j:S[j])
    u=U[:,i]; v=V[i,:]
    c=sum(u[r]*B[r,k]*v[k] for r in range(n) for k in range(n))
    Bn=mpm.mnorm(B,'F')
    # largest eigenvalue of A^{-1}B via generalized: eig of A^{-1}B
    DO_EIG = n <= 30
    M = mpm.inverse(A)*B if DO_EIG else None
    ev = mpm.eig(M,left=False,right=False) if DO_EIG else [mpm.mpf(1)]
    lmax=max(abs(e) for e in ev)
    print((n,b,q,g,h),'-log smin/n=%.4f  log|u^T B v|/n=%.4f  log||B||_F/n=%.4f  log|lambda_max(A^-1 B)|/n=%.4f  [%.1fs]'%(
        -float(mpm.log(S[i]))/n, float(mpm.log(abs(c)))/n, float(mpm.log(Bn))/n, float(mpm.log(lmax))/n, time.time()-t),flush=True)
