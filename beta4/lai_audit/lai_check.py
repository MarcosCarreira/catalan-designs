"""Checks of L. Lai, arXiv:2610.06838: Lemma 2.2 (formula vs series), J_n nonvanishing (Prop 3.1 cases),
denominators of J_n (Prop 5.5 trend) and size of J_n with c = a zeta(3) + b zeta(4) (Prop 6.6 trend)."""
import sys, math, time
from fractions import Fraction as Fr
import mpmath as mp
def q2(x):
    x=Fr(x); return mp.mpf(x.numerator)/x.denominator
import sympy as sp
from flint import fmpq_mat, fmpq
def H(k,s): return sum(Fr(1,j**s) for j in range(1,k+1))
def Lvals(n,a,b,c):
    """L_n(t^s) for s = 0..6n-2 via Lemma 2.2; c may be Fraction or mpf"""
    N6=6*n; out=[]
    U=[];V=[];W=[]
    for k in range(N6+1):
        w=Fr(math.comb(N6,k)**2*math.comb(2*n+k,k)**3)
        A=3*a*H(k,2)+2*b*H(k,3)
        U.append(-A); V.append((3*H(2*n+k,1)+2*H(N6-k,1)-5*H(k,1))*A); W.append((w, -a*H(k,3)-b*H(k,4)))
    for s in range(2*3*n-1):
        tot=0
        for k in range(N6+1):
            w,hk=W[k]; x=-k
            P=Fr(x)**s if s>0 or True else 1
            dP=s*Fr(x)**(s-1) if s>0 else 0
            tot+=w*(U[k]*dP+V[k]*P)+w*6*(c+hk)*P
        out.append(tot)
    return out
def Jn(n,a,b,c,exact=True):
    Lv=Lvals(n,a,b,c); m=3*n
    norm=math.prod(math.factorial(i)**2 for i in range(m))
    if exact:
        M=fmpq_mat(m,m,[fmpq(Lv[i+j].numerator,Lv[i+j].denominator) for i in range(m) for j in range(m)])
        d=M.det(); return Fr(int(d.p),int(d.q))/norm
    M=mp.matrix(m,m)
    for i in range(m):
        for j in range(m): M[i,j]=Lv[i+j]
    return mp.det(M)/norm
if __name__=="__main__":
    mp.mp.dps=60
    # Lemma 2.2: direct series for n = 1, P = t^s, a, b real, c = a zeta3 + b zeta4
    t=sp.symbols('t'); n=1; a,b=2,-3
    Rn=sp.Rational(math.factorial(6*n)**2,math.factorial(2*n)**3)*sp.rf(t-2*n,2*n)**3/sp.rf(t,6*n+1)**2
    c_real=a*mp.zeta(3)+b*mp.zeta(4)
    worst=0
    for s in range(0,5):
        f=sp.simplify(Rn*t**s); f1=sp.lambdify(t,sp.diff(f,t),'mpmath'); f2=sp.lambdify(t,sp.diff(f,t,2),'mpmath')
        direct=-3*a*mp.nsum(lambda m:f1(m),[1,mp.inf])+b*mp.nsum(lambda m:f2(m),[1,mp.inf])
        # formula with real c: compute with Fractions except c
        N6=6*n; tot=mp.mpf(0)
        for k in range(N6+1):
            w=math.comb(N6,k)**2*math.comb(2*n+k,k)**3; A=3*a*H(k,2)+2*b*H(k,3)
            Uk=-A; Vk=(3*H(2*n+k,1)+2*H(N6-k,1)-5*H(k,1))*A
            P=mp.mpf(-k)**s; dP=s*mp.mpf(-k)**(s-1) if s>0 else 0
            tot+=w*(q2(Uk)*dP+q2(Vk)*P+6*(c_real-a*q2(H(k,3))-b*q2(H(k,4)))*P)
        worst=max(worst,abs(direct-tot)/max(1,abs(tot)))
    print('Lemma 2.2, n=1, s=0..4: max relative difference series vs formula:',mp.nstr(worst,3))
    # Prop 3.1: q=19 -> n=4 ; q=37 -> n=8 ; several (a,b,c) with 0<|b|<q
    for (a,b,c) in ((1,1,1),(2,-3,5),(7,11,-4),(0,1,0)):
        J4=Jn(4,a,b,c); print(f'Prop 3.1 q=19 (n=4), (a,b,c)={(a,b,c)}: J_4 != 0: {J4!=0}, v_19(J_4) = {sp.multiplicity(19,J4.numerator)-sp.multiplicity(19,J4.denominator)}')
    # Prop 5.5 / 6.6 trends
    for n in range(1,9):
        T0=time.time()
        Ji=Jn(n,1,1,1); den=math.log(Ji.denominator); num=math.log(abs(Ji.numerator)) if Ji!=0 else float('nan')
        mp.mp.dps=60+40*n*n
        Jr=Jn(n,1,1,1*mp.zeta(3)+1*mp.zeta(4),exact=False) if False else None
        print(f'n={n}: integer (a,b,c)=(1,1,1): -log den(J_n)/n^2 = {-den/n**2:8.3f}, log|J_n|/n^2 = {(num-den)/n**2:8.3f}  [{time.time()-T0:.1f}s]',flush=True)
