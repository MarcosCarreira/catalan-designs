import sys, math
exec(open('sun_check.py').read().split('for B in map')[0])
import itertools
def parts(B):
    S=max(1,B//20); N=2*B+S+3
    Ssum=[fmpq(0)]
    for k in range(N+S+5): Ssum.append(Ssum[-1]+fmpq((-1)**k,(2*k+1)**2))
    u=lambda m:(fmpq_poly([-Ssum[m],1])*((-1)**m))/(2*m+1)
    Pi=[math.prod((2*(h+i)+1)**2 for h in range(1,B+1)) for i in range(N+S+3)]
    R=[[sum((u(i+j)*((-1)**i*math.comb(a+2*B,i)*Pi[i]) for i in range(a+2*B+1)),fmpq_poly([0])) for j in range(1,S+1)] for a in range(S+3)]
    A=tuple(range(S))
    det=fmpq_poly([0])
    for perm in itertools.permutations(range(S)):
        sgn=1
        for x in range(S):
            for y in range(x+1,S):
                if perm[x]>perm[y]: sgn=-sgn
        t=fmpq_poly([sgn])
        for r_,c_ in zip(A,perm): t=t*R[r_][c_]
        det+=t
    mp.mp.dps=4000; Gv=+mp.catalan
    dv=mp.mpf(0)
    for k in range(det.degree(),-1,-1): dv=dv*Gv+mp.mpf(int(det[k].p))/int(det[k].q)
    lF=sum(mp.log(mp.factorial(r)) for r in range(2*B)); lPi=sum(mp.log(Pi[i]) for i in range(N))
    Dd=1
    for k in range(det.degree()+1): Dd=Dd*int(det[k].q)//math.gcd(Dd,int(det[k].q))
    print(f'B={B} S={S}: log|det R(G)|={float(mp.log(abs(dv))):.1f}  log den(det R)={math.log(Dd):.1f}  log F_B={float(lF):.1f}  log prod Pi={float(lPi):.1f}  (per B^2: {float(mp.log(abs(dv)))/B**2:.3f}, {math.log(Dd)/B**2:.3f}, {float(lF)/B**2:.3f}, {float(lPi)/B**2:.3f})')
    # also a single entry R[0][0]: real value vs crude size
    e=R[0][0]; ev=mp.mpf(int(e[0].p))/int(e[0].q)+Gv*mp.mpf(int(e[1].p))/int(e[1].q)
    print(f'   R[0][1] at G: log|.| = {float(mp.log(abs(ev))):.1f};  log max_i C(2B,i) Pi_i |u_i| ~ {float(max(mp.log(math.comb(2*B,i)*Pi[i]) for i in range(2*B+1))-2*mp.log(2*B)):.1f}')
for B in map(int,sys.argv[1].split(',')): parts(B)
