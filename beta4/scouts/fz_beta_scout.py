"""fz_beta_scout.py (scout) -- Fauzan-type Hankel determinant for beta(s), s even, chi_-4 analogue.
Functional: mu(t^e) = (2e+s)! (2^(2e+2)-1) |B_(2e+2)| / ((s-1)! (2e+1)! (2e+2)) (integral against
w(y) = (-1)^(s-1) y^s g^(s-1)(y)/(s-1)!, g = 1/sinh(pi y), alternating Abel-Plana formula);
mu_X(1/(t+(j+1/2)^2)) = 2 (2j+1)^(s-1) (-1)^j (X - sum_(l<j) (-1)^l (2l+1)^-s) - 1/(2j+1), X = beta(s).
Rows D_N(t)^r t^(i+j)/D_K(t), D_m(t) = prod_(j<m) (t+(j+1/2)^2). Q = primitive integer polynomial of det.
usage: python3 fz_beta_scout.py K N h r s [check]"""
import sys, math, time
from flint import fmpq, fmpq_poly, fmpz_mat
import mpmath as mp
K,N,h,r,S=(int(x) for x in sys.argv[1:6]); T0=time.time()
def bern(m): A,B=mp.bernfrac(m); return fmpq(int(A),int(B))
def mu_t(e): return fmpq(math.factorial(2*e+S)*(2**(2*e+2)-1), math.factorial(S-1)*math.factorial(2*e+1)*(2*e+2))*abs(bern(2*e+2))
t=fmpq_poly([0,1])
P=lambda j: t+fmpq((2*j+1)**2,4)
DN=fmpq_poly([1]);  [DN:=DN*P(j) for j in range(N)]
Dt=fmpq_poly([1]);  [Dt:=Dt*P(j) for j in range(N,K)]
Dtp=Dt.derivative()
Bp=[fmpq(0)]
for l in range(K): Bp.append(Bp[-1]+fmpq((-1)**l,(2*l+1)**S))
DNr=DN**(r-1)          # D_N^r / D_K = D_N^(r-1) / D_tail
maxdeg=(r-1)*N+2*h
mut=[mu_t(e) for e in range(maxdeg+1)]
a=[];b=[]
for s in range(2*h-1):
    q,rem=divmod(DNr*t**s,Dt)
    A=sum((q[e]*mut[e] for e in range(q.degree()+1)),fmpq(0)); B=fmpq(0)
    for j in range(N,K):
        c=rem(-fmpq((2*j+1)**2,4))/Dtp(-fmpq((2*j+1)**2,4))
        k=2*(2*j+1)**(S-1)*(-1)**j
        B+=c*k; A+=c*(-k*Bp[j]-fmpq(1,2*j+1))
    a.append(A);b.append(B)
L=1
for q in a+b: L=L*int(q.q)//math.gcd(L,int(q.q))
ai=[int((q*L).p) for q in a]; bi=[int((q*L).p) for q in b]
xs=list(range(h+1)); ys=[fmpq(int(fmpz_mat(h,h,[ai[i+j]+x*bi[i+j] for i in range(h) for j in range(h)]).det())) for x in xs]
coef=ys[:]
for k in range(1,h+1):
    for i in range(h,k-1,-1): coef[i]=(coef[i]-coef[i-1])/(xs[i]-xs[i-k])
Pp=fmpq_poly([coef[h]])
for k in range(h-1,-1,-1): Pp=Pp*(t-xs[k])+coef[k]
assert Pp.degree()==h, Pp.degree()
g=0;D=1
for k in range(h+1): g=math.gcd(g,int(Pp[k].p)); D=D*int(Pp[k].q)//math.gcd(D,int(Pp[k].q))
Q=Pp*fmpq(D,g)
bits=max(abs(int(Q[k].p)).bit_length() for k in range(h+1))
vals=[]
for ex in (300,600):
    mp.mp.dps=int(0.302*bits)+ex
    b4=(mp.zeta(S,mp.mpf(1)/4)-mp.zeta(S,mp.mpf(3)/4))/4**S
    v=mp.mpf(0)
    for k in range(h,-1,-1): v=v*b4+int(Q[k].p)
    vals.append(v)
v=vals[1]; ok=abs(vals[0]-v)/abs(v)<mp.mpf(10)**-30
lq=mp.log(abs(v))
extra=''
if 'check' in sys.argv:
    mp.mp.dps=int(0.302*bits)+300
    b4=(mp.zeta(S,mp.mpf(1)/4)-mp.zeta(S,mp.mpf(3)/4))/4**S
    G=mp.matrix(h,h)
    for i in range(h):
        for j in range(h): G[i,j]=mp.mpf(int(a[i+j].p))/int(a[i+j].q)+b4*mp.mpf(int(b[i+j].p))/int(b[i+j].q)
    dd=mp.det(G); extra=f' det(G(beta4))>0:{dd>0}'
print(f"s={S} K={K} N={N} h={h} r={r}: Q(beta_s) {'>0' if v>0 else '<0'} stable={ok} log|Q|={float(lq):.2f} per K^2={float(lq)/K**2:+.5f} log b-bound={float(-lq/h):+.4f} log maxcoef={bits*0.693:.0f}{extra} [{time.time()-T0:.1f}s]",flush=True)
