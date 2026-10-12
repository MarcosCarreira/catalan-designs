from fractions import Fraction as F
import math, sys
exec(open('consts.py').read().split('# integral of R')[0])
def brk(lo,hi):
    S={F(lo),F(hi)}
    for d in (F(1),F(1,2),1/(2*H),1/(2*lam),1/al,1/(2*al),1/(1-al),1/(1+al),1/(2*(1-al)),1/(2*(1+al))):
        k=math.ceil(F(lo)/d)
        while k*d<=hi: S.add(k*d); k+=1
    return sorted(S)
def piece_int(a,b,depth=0):
    # exact integral of R(x)/x^3 on [a,b] assuming R linear there; verified by collinearity, else bisect
    e=(b-a)/1000
    x0,x1,x2=a+e,(a+b)/2,b-e
    r0,r1,r2=R(x0),R(x1),R(x2)
    if (r1-r0)*(x2-x0)!=(r2-r0)*(x1-x0):
        if depth>30: raise Exception('no linear piece')
        # kink: find by solving s=n+ not needed; bisect
        m=(a+b)/2; return piece_int(a,m,depth+1)+piece_int(m,b,depth+1)
    beta=(r2-r0)/(x2-x0); alpha=r0-beta*x0
    return (-alpha/(2*b*b)-beta/b)-(-alpha/(2*a*a)-beta/a)
def total(lo,hi):
    P=brk(lo,hi); Q=[]
    for a,b in zip(P,P[1:]):
        m=(a+b)/2; T=fl(2*H*m); q=fl(2*m); xs=F(T-q,2)/(H-1)
        Q.append(a)
        if a<xs<b: Q.append(xs)
    Q.append(P[-1])
    return sum(piece_int(a,b) for a,b in zip(Q,Q[1:]))
lo,hi=int(sys.argv[1]),int(sys.argv[2])
v=total(lo,hi)
print(f'int_{lo}^{hi} R/x^3 =',float(v))
if (lo,hi)==(3,20):
    ref=F(322437603634266857629,7535670527041937280000); print('paper:',float(ref),' equal exactly:',v==ref, ' diff', float(v-ref))
