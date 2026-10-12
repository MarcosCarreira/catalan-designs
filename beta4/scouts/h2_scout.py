"""H2 scout: determinants with entries in Q + Q beta(4) from half-integer arrays (even indices).
Entry(r,k) = sum_m sum_{i,j} a^{(m)}_{r,i} c_{k,j} Z'_m(i,j), Z'_m(i,j) = sum_u (-1)^u/((i+u+1/2)^m (j+u+1/2)).
a^{(3)}_r = coefficients of t^{2(alpha+r)} (1-t^2)^h ; columns c_k = coefficients of s^{2(b+k)} (1-s^2)^q.
a^{(2)}, a^{(1)} on the column support J chosen so that G and pi^3 cancel; pi cancels by parity.
Measures: true size log|Delta|/n^2, and the denominator of Delta with beta(4) replaced by 1/3."""
import sys, math, time
from fractions import Fraction as F
from math import comb
import mpmath as mp
T0=time.time()
def Rm(j,m):  # sum_{v<j} (-1)^v/(v+1/2)^m
    return sum(F((-1)**v)/F(2*v+1,2)**m for v in range(j))
cacheR={}
def R(j,m):
    if (j,m) not in cacheR: cacheR[(j,m)]=Rm(j,m)
    return cacheR[(j,m)]
# T_m(j) = (-1)^j (2^m beta(m) - R(j,m)); for even j: 2^m beta(m) - R(j,m)
def Zrat(m,i,j):
    """rational part of Z'_m(i,j) and the coefficients of beta(1..4) (i, j even)."""
    # returns dict {0: rational, 1:..,2:..,3:..,4:..}
    out={0:F(0),1:F(0),2:F(0),3:F(0),4:F(0)}
    def addT(mm,idx,coef):   # coef * T_mm(idx), idx even
        out[mm]+=coef*2**mm; out[0]-=coef*R(idx,mm)
    if i==j:
        addT(m+1,i,F(1)); return out
    d=i-j
    if m==1: addT(1,j,F(1,d)); addT(1,i,F(-1,d))
    if m==2: addT(1,j,F(1,d*d)); addT(1,i,F(-1,d*d)); addT(2,i,F(-1,d))
    if m==3: addT(1,j,F(1,d**3)); addT(1,i,F(-1,d**3)); addT(2,i,F(-1,d*d)); addT(3,i,F(-1,d))
    return out
def build(n,alpha,h,b,q,W=4):
    rows=[]; 
    for r in range(n):
        a3={2*(alpha+r+l):F((-1)**l*comb(h,l)) for l in range(h+1)}
        rows.append(a3)
    cols=[{2*(b+k+l):F((-1)**l*comb(q,l)) for l in range(q+1)} for k in range(n)]
    J=sorted(set(j for c in cols for j in c))
    Frat=[[F(0)]*n for _ in range(n)]; Fb4=[[F(0)]*n for _ in range(n)]; resid=0
    for r,a3 in enumerate(rows):
        # a2_j = sum_{i != j} a3_i/(i-j) ; a1_j = sum_{i != j} a2_i/(i-j) + sum_{i!=j} a3_i/(i-j)^2   (j in J)
        if W==4:
            a2={j: sum((v/(i-j) for i,v in a3.items() if i!=j),F(0)) for j in J}
            a1={j: sum((v/(i-j) for i,v in a2.items() if i!=j),F(0))+sum((v/F(i-j)**2 for i,v in a3.items() if i!=j),F(0)) for j in J}
            arrs=((3,a3),(2,a2),(1,a1))
        else:
            arrs=((1,a3),)
        for k,c in enumerate(cols):
            tot={0:F(0),1:F(0),2:F(0),3:F(0),4:F(0)}
            for j,cv in c.items():
                for m,arr in arrs:
                    for i,av in arr.items():
                        if av==0: continue
                        z=Zrat(m,i,j)
                        for key in tot: tot[key]+=av*cv*z[key]
            tgt=4 if W==4 else 2
            if any(tot[x]!=0 for x in (1,2,3,4) if x!=tgt): resid+=1
            Frat[r][k]=tot[0]; Fb4[r][k]=tot[tgt]
    return Frat,Fb4,resid
def det_frac(M):
    n=len(M); A=[[x for x in row] for row in M]; d=F(1)
    for c in range(n):
        p=next((r for r in range(c,n) if A[r][c]!=0),None)
        if p is None: return F(0)
        if p!=c: A[c],A[p]=A[p],A[c]; d=-d
        d*=A[c][c]; inv=1/A[c][c]
        for r in range(c+1,n):
            f=A[r][c]*inv
            if f: A[r]=[x-f*y for x,y in zip(A[r],A[c])]
    return d
def run(n,ra,rh,rb,rq,W=4):
    alpha,h,b,q=round(ra*n),round(rh*n),round(rb*n),round(rq*n)
    Frat,Fb4,resid=build(n,alpha,h,b,q,W)
    mp.mp.dps=60+4*n*n
    b4=(mp.zeta(4,mp.mpf(1)/4)-mp.zeta(4,mp.mpf(3)/4))/256 if W==4 else mp.catalan
    M=mp.matrix(n,n)
    for r in range(n):
        for k in range(n):
            x=Frat[r][k]; y=Fb4[r][k]
            M[r,k]=mp.mpf(x.numerator)/x.denominator+b4*mp.mpf(y.numerator)/y.denominator
    dtrue=mp.det(M)
    Dh=det_frac([[Frat[r][k]+F(1,3)*Fb4[r][k] for k in range(n)] for r in range(n)])
    Lt=mp.log(abs(dtrue))/n**2 if dtrue!=0 else -mp.inf
    Ld=-mp.log(Dh.denominator)/n**2 if Dh!=0 else None
    nz=sum(1 for row in Fb4 for x in row if x!=0)
    return alpha,h,b,q,resid,nz,Lt,Ld,Dh!=0
if __name__=="__main__":
    ra,rh,rb,rq=(float(x) for x in sys.argv[1:5]); ns=[int(x) for x in sys.argv[5].split(',')]; W=int(sys.argv[6]) if len(sys.argv)>6 else 4
    for n in ns:
        alpha,h,b,q,resid,nz,Lt,Ld,nonzero=run(n,ra,rh,rb,rq,W)
        print(f"n={n} (alpha,h,b,q)=({alpha},{h},{b},{q}): unwanted-constant residues {resid}, beta4 entries {nz}; "
              f"log|Delta|/n^2 = {mp.nstr(Lt,6)}, -log den(Delta_hyp)/n^2 = {mp.nstr(Ld,6) if Ld is not None else 'det 0'}"
              f"  gap(true - low) = {mp.nstr(Lt-Ld,5) if Ld is not None else '-'}  [{time.time()-T0:.0f}s]",flush=True)
