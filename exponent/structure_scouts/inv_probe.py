"""inv_probe.py -- smallest singular value of the Catalan matrix A = Phi0 + G Phi1 (no G-columns separated), via the
exact-ball inverse (python-flint arb): prints log max|A^-1|/n and log||A^-1||_F/n, which bracket -log(sigma_min)/n
within log(n)/n. Settings below; no arguments. Log: inv_probe_log_<time>.txt next to the script.
Sandbox values (preprint ratio 48:7:4:4:2): n = 24: 4.1025, 48: 4.2572, 72: 4.3409."""
import sys, os, time, math
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import flint, full_cat
DESIGNS = [(96, 14, 8, 8, 4), (120, 18, 10, 10, 5), (144, 21, 12, 12, 6)]   # preprint ratio, rounded
STAMP = time.strftime('%Y%m%d_%H%M%S'); LOG = os.path.join(HERE, 'inv_probe_log_' + STAMP + '.txt')
def say(m):
    print(m, flush=True)
    with open(LOG, 'a') as fh: fh.write(m + '\n')
for (n,b,q,g,h) in DESIGNS:
    t=time.time(); Frat,FG=full_cat.build(n,b,q,g,h)
    flint.ctx.prec=int(3.4*n*n)+2000
    G=flint.arb.const_catalan()
    A=flint.arb_mat(n,n)
    for r in range(n):
        for k in range(n):
            x,y=Frat[r][k],FG[r][k]
            A[r,k]=flint.arb(flint.fmpq(x.numerator,x.denominator))+G*flint.arb(flint.fmpq(y.numerator,y.denominator))
    Ai=A.inv()
    mx=max(abs(Ai[i,j]).mid() for i in range(n) for j in range(n))
    lm=float(mx.log()) if hasattr(mx,'log') else math.log(float(mx))
    # frobenius norm of inverse
    fro=sum((Ai[i,j]**2) for i in range(n) for j in range(n)).sqrt()
    lf=float(fro.mid().log())
    say(str((n,b,q,g,h))+'  log max|Ainv|/n=%.4f  log||Ainv||_F/n=%.4f (bounds on -log sigma_min/n within log(n)/n=%.3f)  [%.1fs]'%(lm/n,lf/n,math.log(n)/n,time.time()-t))
say('log saved to '+LOG)
