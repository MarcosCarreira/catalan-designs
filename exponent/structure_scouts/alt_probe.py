import math, sys, numpy as np
ARGS=sys.argv[1:]
exec(open('tp_head.py').read())
for des in [tuple(map(int,a.split(','))) for a in ARGS]:
    prob=Problem(des,2,300,200); v,L,Q,w=F(prob,1.0,None,6000)
    Nx=prob.Nx
    mx=(prob.Xw*(w[:Nx]/prob.hx)[:,None]).ravel(); xs=prob.Xg.ravel()
    ms=(prob.Sw*(w[Nx:]/prob.hs)[:,None]).ravel(); ss=prob.Sg.ravel()
    th=np.linspace(0,2*np.pi,721)
    z=np.exp(1j*th)
    ar=np.array([mx@np.log(np.abs(1+xs*zz)) for zz in z]); ac=np.array([ms@np.log(np.abs(1+ss*zz)) for zz in z])
    s=-(L+2*Q)
    print(des,'-dF=%.4f'%s,' row max_z=%.4f (z=%s)'%(ar.max(),np.round(z[ar.argmax()],2)),' col max_z=%.4f'%ac.max(),
          ' sum=%.4f  -dF+row+col=%.4f  -dF+log2=%.4f'%(ar.max()+ac.max(), s+ar.max()+ac.max(), s+math.log(2)),
          ' mass(mu)=%.3f mass(nu)=%.3f'%(mx.sum(),ms.sum()),flush=True)
