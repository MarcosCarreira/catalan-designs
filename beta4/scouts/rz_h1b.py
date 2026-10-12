exec(open('rz_h1.py').read().split('ETA1=')[0])
m.mp.dps=200
b2=m.catalan; b4=(m.zeta(4,m.mpf(1)/4)-m.zeta(4,m.mpf(3)/4))/256
import math
def lcmden(*xs):
    D=1
    for x in xs: D=D*x.denominator//math.gcd(D,x.denominator)
    return D
for eta in ((3,1,1,1,1,1,1),(5,2,2,2,2,2,2),(7,3,3,3,3,3,3),(8,3,3,3,3,3,3)):
    q=len(eta)-1
    if not (sum(eta[1:])*2<=(q-1)*eta[0] and all(0<2*e<eta[0] for e in eta[1:])): print(eta,'violates (3)'); continue
    print(eta)
    for n in (4,8,12,16):
        (c0,c2,c4),_=full_form(eta,n)
        D=lcmden(c0,c2,c4)
        S=m.mpf(c2.numerator)/c2.denominator*b2+m.mpf(c4.numerator)/c4.denominator*b4+m.mpf(c0.numerator)/c0.denominator
        dn=1
        for k in range(1,n+1): dn=dn*k//math.gcd(dn,k)
        mx=max(abs(c0),abs(c2),abs(c4))
        ints=[int(x*D) for x in (c0,c2,c4)]; g=math.gcd(math.gcd(ints[0],ints[1]),ints[2])
        red=m.log(abs(S))+m.log(D)-m.log(g)
        print(f"   n={n}: log|S|/n={m.nstr(m.log(abs(S))/n,5)}  log D/n={m.nstr(m.log(D)/n,5)} (log d_n/n={m.nstr(m.log(dn)/n,4)})  log(content)/n={m.nstr(m.log(g)/n,5)}  log|integral form|/n={m.nstr(red/n,5)}",flush=True)
