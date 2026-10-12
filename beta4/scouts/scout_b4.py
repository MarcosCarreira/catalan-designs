import mpmath as m, time
T0=time.time()
m.mp.dps=110
pi,L2,G=m.pi,m.log(2),m.catalan
z=(1+1j)/2
ImLi3=m.im(m.polylog(3,m.mpc(0.5,0.5))); ImLi4=m.im(m.polylog(4,m.mpc(0.5,0.5)))
beta4=(m.zeta(4,m.mpf(1)/4)-m.zeta(4,m.mpf(3)/4))/256
z3=m.zeta(3)
def S0(r): return 2/m.factorial(r)*m.quad(lambda t: t*(-m.log(m.sin(t)))**r/m.sin(t),[0,pi/4,pi/2])
def S1(r): return 2/m.factorial(r)*m.quad(lambda t: (t/m.sin(t)**2-m.cot(t))*(-m.log(m.sin(t)))**r,[0,pi/4,pi/2])
def P1(r): return 2/m.factorial(r)*m.quad(lambda t: t*(-m.log(m.sin(t)))**r,[0,pi/4,pi/2])
bases={
 1:(['pi','pi*log2','G'],[pi,pi*L2,G]),
 2:(['G','pi^2','pi*log2','log2^2'],[G,pi**2,pi*L2,L2**2]),
 3:(['pi^3','pi*log2^2','G*log2','ImLi3','zeta3','pi^2*log2','log2^3'],[pi**3,pi*L2**2,G*L2,ImLi3,z3,pi**2*L2,L2**3]),
 4:(['beta4','ImLi4','ImLi3*log2','G*log2^2','pi^2*G','pi^3*log2','pi*log2^3','pi*zeta3','pi^4','zeta3*log2','pi^2*log2^2','log2^4','Li4half'],
    [beta4,ImLi4,ImLi3*L2,G*L2**2,pi**2*G,pi**3*L2,pi*L2**3,pi*z3,pi**4,z3*L2,pi**2*L2**2,L2**4,m.polylog(4,0.5)]),
}
def ident(name,val,w):
    names,vals=bases[w]
    rel=m.pslq([val,1]+vals,maxcoeff=10**8,maxsteps=10**6)
    if rel is None: print(f"  {name}: no relation (weight {w} basis + 1)"); return
    terms=[f"{-c}/{rel[0]}*{n}" for c,n in zip(rel[2:],names) if c]
    print(f"  {name} = {-rel[1]}/{rel[0]} " + " ".join(terms))
for r in range(0,3):
    a,b,c=S0(r),S1(r),P1(r)
    print(f"r={r}: S0={m.nstr(a,20)} S1={m.nstr(b,20)} P1={m.nstr(c,20)}  [{time.time()-T0:.1f}s]",flush=True)
    ident('S0',a,r+2); ident('S1',b,r+1 if r>0 else 1); ident('P1',c,r+2)
