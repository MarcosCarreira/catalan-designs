import mpmath as m, time
T0=time.time()
m.mp.dps=70
pi,L2,G,z3=m.pi,m.log(2),m.catalan,m.zeta(3)
b4=(m.zeta(4,m.mpf(1)/4)-m.zeta(4,m.mpf(3)/4))/256
I4=m.im(m.polylog(4,m.mpc(0.5,0.5))); I3=m.im(m.polylog(3,m.mpc(0.5,0.5)))
Li4h=m.polylog(4,m.mpf(1)/2)
# A_n = sum 4^k/((2k+1)^n C(2k,k)),  B_n = sum 4^k/((2k+2)^n C(2k,k)),  C_n = sum_{k>=1} 4^k/((2k)^n C(2k,k))
def chi(n,y): return (m.polylog(n,y)-m.polylog(n,-y))/2
def A(n): return m.quad(lambda t: chi(n-1,m.sin(t)),[0,pi/4,pi/2])
def B(n): return m.quad(lambda t: m.polylog(n,m.sin(t)**2)/(2**n*m.sin(t)),[0,pi/4,pi/2])
def C(n): return m.quad(lambda t: m.sin(t)*m.polylog(n,m.sin(t)**2)/2**n,[0,pi/4,pi/2])
names=['1','beta4','ImLi4','ImLi3*log2','G*log2^2','pi^2*G','pi^3*log2','pi*log2^3','pi*zeta3','pi^4','zeta3*log2','pi^2*log2^2','log2^4','Li4half','pi','pi*log2','pi^2','log2','zeta3','pi^3','G','pi^2*log2','pi*log2^2','log2^2','log2^3','ImLi3','G*log2']
vals=[1,b4,I4,I3*L2,G*L2**2,pi**2*G,pi**3*L2,pi*L2**3,pi*z3,pi**4,z3*L2,pi**2*L2**2,L2**4,Li4h,pi,pi*L2,pi**2,L2,z3,pi**3,G,pi**2*L2,pi*L2**2,L2**2,L2**3,I3,G*L2]
def ident(nm,v,sel):
    ns=[names[i] for i in sel]; vs=[vals[i] for i in sel]
    rel=m.pslq([v]+vs,maxcoeff=10**6,maxsteps=10**6)
    if rel is None: print(f"  {nm}: no relation over {ns}",flush=True); return
    print(f"  {nm} = "+" + ".join(f"{m.mpf(-c)/rel[0]}*{n}" for c,n in zip(rel[1:],ns) if c),flush=True)
w4=[0,1,2,3,4,5,6,7,8,9,10,11,12,13]
for nm,f in ():
    v=f(4); print(nm,m.nstr(v,25),f"[{time.time()-T0:.1f}s]",flush=True)
    ident(nm,v,w4)
print("--- even mixed-weight basis")
ev=[0,16,18,21,17,23,24,9,10,11,12,13]
for nm,f,n in (('C2',C,2),('C3',C,3),('C4',C,4),('B3',B,3),('B4',B,4)):
    ident(nm,f(n),ev)
