import math, sys, time, numpy as np
ARGS=sys.argv[1:]
exec(open('tp_head.py').read())
for des in [tuple(map(int,a.split(','))) for a in ARGS]:
    t=time.time(); prob=Problem(des,2,300,200); v,L,Q,w=F(prob,1.0,None,6000)
    s=-(L+2*Q)
    print(des,'F(1) U-conv %.5f'%(v-0.5*math.log(2)),' -dF/dtheta=%.4f'%s,' +log2=%.4f'%(s+math.log(2)),' [%.1fs]'%(time.time()-t),flush=True)
