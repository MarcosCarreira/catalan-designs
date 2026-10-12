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
