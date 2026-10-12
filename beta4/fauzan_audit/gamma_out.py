from fractions import Fraction as F
exec(open('iout.py').read().split('# piecewise linear')[0])
def gamma_out(K,N,p):
    v=K-p*(K//p); u=max(0,N+v-p+1); tp=min(N,v)+u; rp=max(0,K+4*N-2*p+2)
    return (-7*(K-p)+6*tp-1-min(rp,p-1-N+u)) if K<2*p else (-7*(K-p)+3+12*N+5*tp-min(rp,p+u))
from sympy import nextprime
K=320000; N=3*K//40
p=nextprime(K//3)
while p<=K:
    y=F(p,K); pred=K*(R0(y)-d(y))
    print(f'p={p} y={float(y):.3f}: -gamma_out={-gamma_out(K,N,p)}  K(R0-d)={float(pred):.1f}  diff={float(-gamma_out(K,N,p)-pred):.1f}')
    p=nextprime(p+K//9)
