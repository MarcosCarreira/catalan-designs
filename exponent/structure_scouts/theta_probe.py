import math, time, sys, numpy as np
sys.argv=[sys.argv[0]]
import importlib.util
spec=importlib.util.spec_from_file_location('ce','cat_equilibrium.py')
src=open('cat_equilibrium.py').read()
# load only definitions (stop before main run)
cut=src.index('def paper_seq')
import os
ns={'__file__':os.path.abspath('cat_equilibrium.py')}
exec(src[:cut].replace("LOGF = open(LOGPATH, 'w')","LOGF = open('/dev/null','w')"),ns)
Problem=ns['Problem']; ITERS=ns['ITERS']
def F(prob, theta, w0=None, iters=6000):
    M0=prob.M
    prob.M=theta*M0; prob.Lc=abs(np.linalg.eigvalsh(prob.M).min())
    _, w = prob.maximise(0.0, iters, w0)
    prob.M=M0; prob.Lc=abs(np.linalg.eigvalsh(M0).min())
    l=np.concatenate([prob.lx0, prob.ls])
    L=(1+prob.al+prob.gg)*math.log(2)+l@w
    Q=0.5*math.log(2)+0.5*w@M0@w
    return theta*L+theta**2*Q, L, Q, w
for name,des in (('paper',(48,7,4,4,2)),('new',(96,15,6,8,5))):
    t=time.time(); prob=Problem(des,2,300,200)
    vals={}; w=None
    for th in (1.0,0.99,0.98,0.95):
        v,L,Q,w=F(prob,th,w); vals[th]=(v,L,Q)
    v1,L1,Q1=vals[1.0]
    print(name,'F(1) raw=%.5f (U-conv %.5f)'%(v1,v1-0.5*math.log(2)),' L*=%.4f Q*=%.4f'%(L1,Q1))
    print('  envelope F\'(1)=L*+2Q*=%.4f'%(L1+2*Q1), ' finite diff (1,.99): %.4f  (1,.98): %.4f  (1,.95): %.4f'%((v1-vals[0.99][0])/0.01,(v1-vals[0.98][0])/0.02,(v1-vals[0.95][0])/0.05),' [%.1fs]'%(time.time()-t),flush=True)
