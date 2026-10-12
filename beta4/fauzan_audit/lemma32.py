"""Fauzan Lemmas 3.1-3.2: kappa_d valuations; distribution formula, checked p-adically (truncated far-pole series)."""
from fractions import Fraction as F
from sympy import bernoulli, binomial
import math, random
def B(k): return F(bernoulli(k)) if k!=1 else F(-1,2)
def vp(q,p):
    q=F(q)
    if q==0: return 10**9
    a,b=q.numerator,q.denominator; v=0
    while a%p==0: a//=p; v+=1
    while b%p==0: b//=p; v-=1
    return v
Bs=[B(k) for k in range(600)]
# Lemma 3.1 premise: v_p(kappa_d) >= -1 always, kappa_d in Z_p for d <= p+1
bad=[]
for p in (7,11,13,17,19,23,29,31,37,41,43,47):
    for d in range(3,500):
        kd=F(d*(d-1)*(d-2))*Bs[d-3]/24
        if vp(kd,p)<-1 or (d<=p+1 and vp(kd,p)<0): bad.append((p,d,vp(kd,p)))
print('kappa_d checks (p 7..47, d<500): violations',bad[:5],len(bad))
# tau on polynomials: tau(P)=L(P''')/24, L(x^k)=B_k
def tau_poly(c):  # c[k] coefficient of x^k
    return sum(c[d]*F(d*(d-1)*(d-2))*Bs[d-3]/24 for d in range(3,len(c)))
def poly_sub(c,a,p):  # coefficients of P(a+p x)
    out=[F(0)]*len(c)
    for d,cd in enumerate(c):
        if cd==0: continue
        for j in range(d+1): out[j]+=cd*math.comb(d,j)*F(a)**(d-j)*F(p)**j
    return out
random.seed(1); okpoly=True
for p in (7,11,13):
    for trial in range(5):
        c=[F(random.randint(-9,9)) for _ in range(12)]
        lhs=tau_poly(c); rhs=sum(tau_poly(poly_sub(c,a,p)) for a in range(p))/F(p)**4
        okpoly&= lhs==rhs
print('distribution formula on polynomials (exact):',okpoly)
def H5(n): return sum(F(1,v**5) for v in range(1,n+1))
def dd(r): return r if r>=0 else -r-1
KMAX=40
def tau_an_far(s):  # (3.5), truncated
    return -F(1,4)*sum(math.comb(k+3,3)*Bs[k]*F(s)**(-k-4) for k in range(KMAX))
for p in (7,11,13):
    Cp=sum(tau_an_far(F(-a,p)) for a in range(1,p))
    worst=10**9
    for r in list(range(-25,26)):
        s0=None; far=F(0)
        for a in range(p):
            if (r-a)%p==0: s0=(r-a)//p
            else: far+=tau_an_far(F(r-a,p))
        # identity: p^5 H_{d(r)} = H_{d(s0)} - C_p + far  (X-coefficients agree: -p^5 X both sides)
        diff=F(p)**5*H5(dd(r))-(H5(dd(s0))-Cp+far)
        worst=min(worst,vp(diff,p))
    print(f'p={p}: distribution formula at poles r=-25..25: min v_p(difference) = {worst} (truncation level about {KMAX+3})')
# C_p in p^5 Z_p (Section 3.2)
for p in (7,11,13,17,19,23,29,31,37,41,43,47):
    Cp=sum(tau_an_far(F(-a,p)) for a in range(1,p))
    print(f'C_p check p={p}: v_p(C_p) = {vp(Cp,p)}')
