import sys, time, measure_scout as m, mpmath as mpm
def lam(d):
    Z=m.delta_poly(*d)
    if isinstance(Z,tuple): Z=Z[0]
    n=d[0]; digits=len(str(max(abs(int(c)) for c in Z.coeffs())))
    c=m.taylor_at_G(Z, digits+n*n//4+60)
    l=[mpm.log(abs(x)) for x in c]
    s=[float(l[k+1]-l[k]) for k in range(len(c)-1)]
    K=len(s)
    # least squares slope of s_k vs k over middle 80%
    a,b=int(0.1*K),int(0.9*K)
    xs=list(range(a,b)); ys=s[a:b]
    mx=sum(xs)/len(xs); my=sum(ys)/len(ys)
    sl=sum((x-mx)*(y-my) for x,y in zip(xs,ys))/sum((x-mx)**2 for x in xs)
    return -sl, s[0]/n, len(c)-1
for d in [tuple(map(int,a.split(','))) for a in sys.argv[1:]]:
    t=time.time(); L,r,deg=lam(d)
    print(d,'deg',deg,'lambda(step)=%.3f'%L,' first/n=%.3f'%r,' [%.1fs]'%(time.time()-t),flush=True)
