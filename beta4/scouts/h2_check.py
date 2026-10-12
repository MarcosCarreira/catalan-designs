exec(open('h2_scout.py').read().split('if __name__')[0])
import mpmath as mp
mp.mp.dps=40
b4=(mp.zeta(4,mp.mpf(1)/4)-mp.zeta(4,mp.mpf(3)/4))/256
n,alpha,h,b,q=3,1,2,0,2
Frat,Fb4,resid=build(n,alpha,h,b,q)
# direct: entry = sum_u (-1)^u sum_{m,i,j} a^m_i c_j / ((i+u+1/2)^m (j+u+1/2))
rows=[{2*(alpha+r+l):(-1)**l*comb(h,l) for l in range(h+1)} for r in range(n)]
cols=[{2*(b+k+l):(-1)**l*comb(q,l) for l in range(q+1)} for k in range(n)]
J=sorted(set(j for c in cols for j in c))
for r in (0,2):
    a3={i:F(v) for i,v in rows[r].items()}
    a2={j: sum((v/(i-j) for i,v in a3.items() if i!=j),F(0)) for j in J}
    a1={j: sum((v/(i-j) for i,v in a2.items() if i!=j),F(0))+sum((v/F(i-j)**2 for i,v in a3.items() if i!=j),F(0)) for j in J}
    for k in (0,2):
        def term(u):
            s=mp.mpf(0)
            for j,cv in cols[k].items():
                for m,arr in ((3,a3),(2,a2),(1,a1)):
                    for i,av in arr.items():
                        s+=mp.mpf(av.numerator)/av.denominator*cv/((i+u+mp.mpf(1)/2)**m*(j+u+mp.mpf(1)/2))
            return (-1)**int(u)*s
        direct=mp.nsum(term,[0,mp.inf])
        dec=mp.mpf(Frat[r][k].numerator)/Frat[r][k].denominator+b4*mp.mpf(Fb4[r][k].numerator)/Fb4[r][k].denominator
        print(r,k,mp.nstr(direct,25),mp.nstr(dec,25),mp.nstr(direct-dec,3))
