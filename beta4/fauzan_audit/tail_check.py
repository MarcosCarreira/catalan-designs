from fractions import Fraction as F
import math, random
exec(open('consts.py').read().split('# integral of R')[0])
def frac(x): return x-math.floor(x)
def Fx(x): f=frac(x); g=frac(al*x); return 4*lam+2*lam*f-12*lam*g
def P(x): f=frac(x); g=frac(al*x); return 74*g*(1-g)-lam*f*(1-f)
def G0(v): return v*(1-v)*(2*v-1)/6
def C(x): return G0(frac(al*x))*74/al-lam*G0(frac(x))
random.seed(3); qmin=10; qmax=-10; dmax=0; cmax=0
for _ in range(3000):
    x=F(random.randint(20*10**6,400*10**6),10**6)
    q=R(x)-x*Fx(x); qmin=min(qmin,q); qmax=max(qmax,q)
    # P' = F + lam (derivative by small difference where smooth)
    e=F(1,10**7); dP=(P(x+e)-P(x-e))/(2*e); dmax=max(dmax,abs(dP-(Fx(x)+lam)))
    cmax=max(cmax,abs(C(x)))
print('Q(x)=R-xF on 3000 random x in [20,400]: min',float(qmin),' max',float(qmax),' (paper: -1/2 <= Q <= 13/8)')
print("P' - (F+lam): max",float(dmax),'(away from breakpoints should be ~0)')
print('max |C| sampled',float(cmax),' (paper bound < 16)')
# mean of P over a period (period of P is lcm of 1 and 1/al = 40/3 -> 40)
S=F(0); n=40*400
for k in range(n): S+=P(F(40*k,n)+F(1,2*n))
print('mean of P (midpoint rule, 16000 pts) =',float(S/n),' paper P-bar = 2923/240 =',2923/240)
print('P(20) =',P(F(20)),' C(20) =',C(F(20)),' P(200)=',P(F(200)),' C(200)=',C(F(200)))
# C' = P - Pbar
x=F(3141,100); e=F(1,10**7); print("C'(x) vs P(x)-Pbar at x=31.41:",float((C(x+e)-C(x-e))/(2*e)),float(P(x)-F(2923,240)))
