from fractions import Fraction as F
import math
al=F(3,40); lam=F(37,40)
def pos(v): return v if v>0 else F(0)
def R0(y):
    if F(1,3)<y<F(1,2): return 8-9*y-8*al-5*min(al,1-2*y)-5*pos(1+al-3*y)
    if F(1,2)<y<1: return 7*(1-y)-6*min(al,1-y)-6*pos(1+al-2*y)+pos(1+4*al-2*y)
    return F(0)
def d(y): return pos(1+4*al-3*y-pos(1+al-3*y)) if F(1,3)<y<F(1,2) else F(0)
def integrand(y): return R0(y)-d(y)-2*lam*math.floor(1/y)+sum(pos(2*lam-j*y) for j in range(1,6))
# piecewise linear: breakpoints
B={F(1,3),F(1,2),F(1),2*lam}
for c in [(1-al)/2,(1+al)/3,1-al,(1+al)/2,(1+4*al)/2,(1+4*al)/3]: B.add(F(c))
for k in range(1,4): B.add(F(1,k))
for j in range(1,6): B.add(2*lam/j)
B=sorted(b for b in B if F(1,3)<=b<=2*lam)
def pint(a,b,depth=0):
    e=(b-a)/1000; x0,x1,x2=a+e,(a+b)/2,b-e; r0,r1,r2=integrand(x0),integrand(x1),integrand(x2)
    if (r1-r0)*(x2-x0)!=(r2-r0)*(x1-x0):
        m=(a+b)/2; return pint(a,m,depth+1)+pint(m,b,depth+1)
    return (b-a)*(integrand(x0)+ (r2-r0)/(x2-x0)*((a+b)/2-x0))
tot=sum(pint(a,b) for a,b in zip(B,B[1:]))
print('Iout =',tot,float(tot),' paper 127751/96000:',tot==F(127751,96000))
Bd=sorted(b for b in B if b<=F(1,2))
dint=sum((b-a)*d((a+b)/2) for a,b in zip(Bd,Bd[1:]))  # d linear on pieces -> midpoint exact
print('int d =',dint,' paper 9/640:',dint==F(9,640))
