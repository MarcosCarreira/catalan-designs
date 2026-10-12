import math, sys
from fractions import Fraction as Fr
from sympy import primerange
exec(open('lai_sum.py').read().split('for n in map')[0])
def frac(x): return x-math.floor(x)
def Fv(rho,m,n,p):
    f=lambda x:x//p
    return m*m+(2*f(6*n)-2*f(6*n-rho)+3*f(2*n+rho)-3*f(2*n)-5)*m+3*max(0,m-2*f(6*n-rho))+(1 if (m==1 and p>max(6*n-rho,2*n+rho)) else 0)
# (a) Table 1 integral
T=[(1,Fr(8,7),-8,5),(Fr(8,7),Fr(6,5),-32,26),(Fr(6,5),Fr(4,3),-8,1),(Fr(4,3),Fr(3,2),-32,19),(Fr(3,2),Fr(8,5),-2,-5),(Fr(8,5),2,-26,10),(2,Fr(8,3),-1,-1),(Fr(8,3),3,-25,8)]
I=sum((b-a)*c0+(b*b-a*a)/2*c1 for a,b,c0,c1 in T)
print('Table 1 integral =',I,float(I),' paper -657/70:',I==Fr(-657,70))
case=[(3,Fr(7,2),5,-4),(Fr(7,2),4,-2,-2),(4,Fr(9,2),-18,2),(Fr(9,2),6,-54,10)]
print('Lemma 5.4 integrals:',[ (b-a)*c0+(b*b-a*a)/2*c1 for a,b,c0,c1 in case], ' sum',sum((b-a)*c0+(b*b-a*a)/2*c1 for a,b,c0,c1 in case))
print('total -5.689676 -657/70 -63/4 =',-5.689676-657/70-63/4)
def tab(n,p):
    x=Fr(p,n)
    rows=[(1,Fr(8,7),lambda:-8*n+5*p-5),(Fr(8,7),Fr(6,5),lambda:-32*n+26*p-8),(Fr(6,5),Fr(4,3),lambda:-8*n+p-3),(Fr(4,3),Fr(3,2),lambda:-32*n+19*p-6),
          (Fr(3,2),Fr(8,5),lambda:-2*n-5*p-1),(Fr(8,5),2,lambda:-26*n+10*p-4),(2,Fr(8,3),lambda:-n-p-2),(Fr(8,3),3,lambda:-25*n+8*p-5)]
    for a,b,g in rows:
        if a<x<=b: return g()
def sec5(n,p):
    if p*p<=8*n: return None
    if p<=n: return 9*n*frac(2*n/p)-6*n*frac(6*n/p)-6*n-Fr(29,16)*p-3*math.sqrt(n)
    if p<=3*n: return tab(n,p)
    x=p/n
    if x<=3.5: return 5*n-4*p
    if x<=4: return -2*n-2*p
    if x<=4.5: return -18*n+2*p-2
    return -54*n+10*p-10
# (b) finite checks of Lemma 5.4 and (d) per-prime comparison sec5 <= Lemma 4.4 bound
for n in map(int,sys.argv[1].split(',')):
    bad=[]; badfin=[]
    for p in primerange(2,6*n+1):
        s=sec5(n,p)
        if s is None: continue
        b=bound44(n,p)
        if s>b+1e-9: bad.append((p,float(s),b))
        if p>3*n:
            x=p/n
            for rho in range(p):
                M=2+2*((6*n-rho)//p)
                for m in range(M+1):
                    Fm=Fv(rho,m,n,p)
                    if x<=3.5: lb=-m-4*(rho<p-2*n)
                    elif x<=4: lb=-2*m-2*(rho<p-2*n)
                    elif x<=4.5: lb=-2*m-2*(rho<=6*n-p)
                    else: lb=2*m-10*(rho<=6*n-p)
                    if Fm<lb: badfin.append((p,rho,m,Fm,lb))
    print(f'n={n}: primes where the Section-5 bound exceeds the Lemma 4.4 bound: {bad[:5]} ({len(bad)}); finite-check failures in Lemma 5.4: {badfin[:5]} ({len(badfin)})')
