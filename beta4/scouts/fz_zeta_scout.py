"""fz_zeta_scout.py (scout) -- Fauzan's Hankel determinant for zeta(s), general s and parameters (s = 5 is Fauzan's).
mu(t^e) = (2e+s)! |B_(2e+2)| / ((s-1)! (2e+2)!), mu_X(1/(t+j^2)) = j^(s-1) (X - H_j^(s)) - 1/(s-1) + 1/(2j);
rows D_N(t)^r t^(i+j)/D_K(t), D_m(t) = prod_(j<=m) (t+j^2). Q = primitive integer polynomial of det; value at zeta(s).
optional 6th argument q: only the poles -(1+q j)^2, j = 0..K-1 (q = 1: Fauzan).
usage: python3 fz_zeta_scout.py K N h r s [q] [check]"""
import sys, math, time
from flint import fmpq, fmpq_poly, fmpz_mat
import mpmath as mp
K,N,h,r,S=(int(x) for x in sys.argv[1:6]); T0=time.time()
QS=int(sys.argv[6]) if len(sys.argv)>6 and sys.argv[6].isdigit() else 1   # pole step: poles at -(1+QS*j)^2
M=lambda j: 1+QS*j
def bern(m): A,B=mp.bernfrac(m); return fmpq(int(A),int(B))
def mu_t(e): return fmpq(math.factorial(2*e+S), math.factorial(S-1)*math.factorial(2*e+2))*abs(bern(2*e+2))
t=fmpq_poly([0,1])
P=lambda j: t+M(j)**2
DN=fmpq_poly([1]);  [DN:=DN*P(j) for j in range(N)]
Dt=fmpq_poly([1]);  [Dt:=Dt*P(j) for j in range(N,K)]
Dtp=Dt.derivative()
Bp=[fmpq(0)]
for l in range(M(K-1)): Bp.append(Bp[-1]+fmpq(1,(l+1)**S))
DNr=DN**(r-1)          # D_N^r / D_K = D_N^(r-1) / D_tail
maxdeg=(r-1)*N+2*h
mut=[mu_t(e) for e in range(maxdeg+1)]
a=[];b=[]
for s in range(2*h-1):
    q,rem=divmod(DNr*t**s,Dt)
    A=sum((q[e]*mut[e] for e in range(q.degree()+1)),fmpq(0)); B=fmpq(0)
    for j in range(N,K):
        jj=M(j); c=rem(-jj*jj)/Dtp(-jj*jj)
        k=jj**(S-1)
        B+=c*k; A+=c*(-k*Bp[jj]-fmpq(1,S-1)+fmpq(1,2*jj))
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
assert Pp.degree()>=1, Pp.degree()
g=0;D=1
for k in range(h+1): g=math.gcd(g,int(Pp[k].p)); D=D*int(Pp[k].q)//math.gcd(D,int(Pp[k].q))
Q=Pp*fmpq(D,g)
bits=max(abs(int(Q[k].p)).bit_length() for k in range(h+1))
vals=[]
for ex in (300,600):
    mp.mp.dps=int(0.302*bits)+ex
    b4=mp.zeta(S)
    v=mp.mpf(0)
    for k in range(h,-1,-1): v=v*b4+int(Q[k].p)
    vals.append(v)
v=vals[1]; ok=abs(vals[0]-v)/abs(v)<mp.mpf(10)**-30
lq=mp.log(abs(v))
extra=''
if 'check' in sys.argv:
    mp.mp.dps=int(0.302*bits)+300
    b4=mp.zeta(S)
    G=mp.matrix(h,h)
    for i in range(h):
        for j in range(h): G[i,j]=mp.mpf(int(a[i+j].p))/int(a[i+j].q)+b4*mp.mpf(int(b[i+j].p))/int(b[i+j].q)
    dd=mp.det(G); extra=f' det(G(beta4))>0:{dd>0}'
print(f"step={QS} zeta s={S} deg={Pp.degree()} K={K} N={N} h={h} r={r}: Q(zeta_s) {'>0' if v>0 else '<0'} stable={ok} log|Q|={float(lq):.2f} per K^2={float(lq)/K**2:+.5f} log b-bound={float(-lq/h):+.4f} log maxcoef={bits*0.693:.0f}{extra} [{time.time()-T0:.1f}s]",flush=True)
