import mpmath as m, time
T0=time.time()
m.mp.dps=110
pi,L2,G=m.pi,m.log(2),m.catalan
z=(1+1j)/2
ImLi3=m.im(m.polylog(3,m.mpc(0.5,0.5))); ImLi4=m.im(m.polylog(4,m.mpc(0.5,0.5)))
beta4=(m.zeta(4,m.mpf(1)/4)-m.zeta(4,m.mpf(3)/4))/256
z3=m.zeta(3)
bases={
 1:(['pi','pi*log2','G'],[pi,pi*L2,G]),
 2:(['G','pi^2','pi*log2','log2^2'],[G,pi**2,pi*L2,L2**2]),
 3:(['pi^3','pi*log2^2','G*log2','ImLi3','zeta3','pi^2*log2','log2^3'],[pi**3,pi*L2**2,G*L2,ImLi3,z3,pi**2*L2,L2**3]),
 4:(['beta4','ImLi4','ImLi3*log2','G*log2^2','pi^2*G','pi^3*log2','pi*log2^3','pi*zeta3','pi^4','zeta3*log2','pi^2*log2^2','log2^4','Li4half'],
    [beta4,ImLi4,ImLi3*L2,G*L2**2,pi**2*G,pi**3*L2,pi*L2**3,pi*z3,pi**4,z3*L2,pi**2*L2**2,L2**4,m.polylog(4,0.5)]),
}
import mpmath as m

def ident(name,val,names,vals):
    rel=m.pslq([val,1]+vals,maxcoeff=10**8,maxsteps=10**6)
    if rel is None: print(f"  {name}: no relation"); return
    print(f"  {name} = " + " + ".join(f"({m.mpf(-c)/rel[0]})*{n}" for c,n in zip(rel[1:],['1']+names) if c))
n4,v4=bases[4]
for a in (1,2,3):
    for r in (0,1):
        val=m.quad(lambda t: t**a*(-m.log(m.sin(t)))**r/m.sin(t),[0,pi/4,pi/2])
        w=a+1+r
        if w==4: ident(f"int theta^{a} (-log sin)^{r}/sin", val, n4, v4)
        elif w==3: ident(f"int theta^{a} (-log sin)^{r}/sin", val, *bases[3])
        elif w==2: ident(f"int theta^{a} (-log sin)^{r}/sin", val, *bases[2])
# mixed-weight basis for S1(r=2)
S1=2/2*m.quad(lambda t: (t/m.sin(t)**2-m.cot(t))*(-m.log(m.sin(t)))**2,[0,pi/4,pi/2])
ident("S1(r=2)",S1,['pi^2','log2','pi^2*log2','zeta3','log2^2','log2^3','pi*G','G'],[pi**2,L2,pi**2*L2,z3,L2**2,L2**3,pi*G,G])
