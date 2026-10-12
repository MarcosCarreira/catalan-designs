import math
from flint import fmpq, fmpq_mat, fmpq_poly
B=20; S=1; N=2*B+S+3; D=2*B
Ssum=[fmpq(0)]
for k in range(N+10): Ssum.append(Ssum[-1]+fmpq((-1)**k,(2*k+1)**2))
Pi=[math.prod((2*(h+i)+1)**2 for h in range(1,B+1)) for i in range(N+5)]
def u(m,G): return ((-1)**m)*(G-Ssum[m])/(2*m+1)
FB=1
for r in range(2*B): FB*=math.factorial(r)
for A in ([0],[1],[2],[3]):
    Ac=[c for c in range(S+3) if c not in A]
    vals=[]
    for G in (fmpq(0),fmpq(1)):
        rows=[]
        for i in range(N):
            row=[fmpq(i**r) for r in range(D)]+[Pi[i]*u(i+1,G)]+[fmpq(math.comb(i,D+c)) for c in Ac]
            rows.append(row)
        M=fmpq_mat(N,N,[x for r in rows for x in r]); vals.append(M.det())
    # R[a][1] at G=0,1
    Rv=[sum(((-1)**i*math.comb(A[0]+2*B,i)*Pi[i]*u(i+1,G) for i in range(A[0]+2*B+1)),fmpq(0)) for G in (fmpq(0),fmpq(1))]
    print('A=',A,' det(A~)/(F_B det R) at G=0,1:', vals[0]/(FB*Rv[0]), vals[1]/(FB*Rv[1]))
