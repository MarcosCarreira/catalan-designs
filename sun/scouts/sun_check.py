"""Z.-W. Sun, arXiv:2609.04176: the scalar qhat_B = +-F_B det R[A,J] / prod_{i<N} Pi_i as a polynomial in G
(T_m = (-1)^m (G - S_{m-1}), u_m = T_m/(2m+1)); measure log(den of qhat as polynomial) + log|qhat(G)| against the claim
<= -delta0 B^2 + o(B^2), delta0 > 0.00966."""
import sys, math, time, itertools
from flint import fmpq, fmpq_poly, fmpq_mat
import mpmath as mp
def run(B, S=None):
    S = S if S is not None else max(1, B//20)
    N = 2*B+S+3
    X = fmpq_poly([0,1])
    Ssum=[fmpq(0)]
    for k in range(N+S+5): Ssum.append(Ssum[-1]+fmpq((-1)**k,(2*k+1)**2))
    def u(m):   # T_m/(2m+1), T_m = (-1)^m (G - S_{m-1}); S_{m-1} = Ssum[m]
        return (fmpq_poly([-Ssum[m],1])*((-1)**m))/(2*m+1)
    Pi=[]
    for i in range(N+S+3):
        Pi.append(math.prod((2*(h+i)+1)**2 for h in range(1,B+1)))
    R=[[None]*S for _ in range(S+3)]
    for a in range(S+3):
        n=a+2*B
        for j in range(1,S+1):
            acc=fmpq_poly([0])
            for i in range(n+1):
                acc+=u(i+j)*((-1)**i*math.comb(n,i)*Pi[i])
            R[a][j-1]=acc
    FB=1
    for r in range(2*B): FB*=math.factorial(r)
    PiN=math.prod(Pi[i] for i in range(N))
    best=None; out=[]
    for A in itertools.combinations(range(S+3),S):
        # determinant of the S x S polynomial matrix by Leibniz (S small) 
        det=fmpq_poly([0])
        for perm in itertools.permutations(range(S)):
            sgn=1
            for x in range(S):
                for y in range(x+1,S):
                    if perm[x]>perm[y]: sgn=-sgn
            t=fmpq_poly([sgn])
            for r_,c_ in zip(A,perm): t=t*R[r_][c_]
            det+=t
        if det==0: continue
        q=det*fmpq(FB,PiN)
        D=1
        for k in range(q.degree()+1):
            if q[k]!=0: D=D*int(q[k].q)//math.gcd(D,int(q[k].q))
        bits=max(abs(int((q[k]*D).p)).bit_length() for k in range(q.degree()+1))
        mp.mp.dps=int(0.302*bits*1.3)+200; Gv=+mp.catalan
        v=mp.mpf(0)
        for k in range(q.degree(),-1,-1): v=v*Gv+mp.mpf(int((q[k]*D).p))
        mp.mp.dps=int(0.302*bits*1.3)+1200; Gv=+mp.catalan
        v2=mp.mpf(0)
        for k in range(q.degree(),-1,-1): v2=v2*Gv+mp.mpf(int((q[k]*D).p))
        assert abs(v-v2)<=abs(v2)*mp.mpf(10)**-20, 'precision'
        L=float(mp.log(abs(v2)))   # = log D + log|qhat(G)|
        out.append((A,L))
        if best is None or L<best[1]: best=(A,L)
    return S,N,out,best
for B in map(int,sys.argv[1].split(',')):
    T=time.time(); S,N,out,best=run(B)
    print(f'B={B} S={S} N={N}: {len(out)} admissible A; best log(H)+log|qhat| = {best[1]:.2f} (per B^2 {best[1]/B**2:+.4f}); worst per B^2 {max(L for _,L in out)/B**2:+.4f}; claim: <= -0.00966 B^2 + o(B^2) = {-0.00966*B*B:.1f} [{time.time()-T:.1f}s]',flush=True)
