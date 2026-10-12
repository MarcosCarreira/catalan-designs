import mpmath as m
from fractions import Fraction as F
from math import comb
m.mp.dps=80
pi,L2,G,z3=m.pi,m.log(2),m.catalan,m.zeta(3)
b4=(m.zeta(4,m.mpf(1)/4)-m.zeta(4,m.mpf(3)/4))/256
I4=m.im(m.polylog(4,m.mpc(0.5,0.5))); I3=m.im(m.polylog(3,m.mpc(0.5,0.5)))
# J_n = int_0^{pi/2} theta^2 sin^n = r_n + pi p_n (n odd): n^2 J_n = n(n-1) J_{n-2} + pi - 2 W_n
W={1:F(1)}; p={1:F(1)}; r={1:F(-2)}
for n in range(3,40,2):
    W[n]=W[n-2]*F(n-1,n); p[n]=F(n-1,n)*p[n-2]+F(1,n*n); r[n]=F(n-1,n)*r[n-2]-2*W[n]/(n*n)
ok=all(abs(m.quad(lambda t: t**2*m.sin(t)**n,[0,pi/2]) - (m.mpf(r[n].numerator)/r[n].denominator + pi*m.mpf(p[n].numerator)/p[n].denominator)) < m.mpf(10)**-60 for n in (1,3,5,11,21))
print("recurrence check (Q + Q pi):", ok)
# closed form of the pi-part: p_{2k+1} = W_{2k+1} * sum_{j<=k} 1/((2j+1)^2 W_{2j+1})
for k in range(6):
    n=2*k+1; s=sum(F(1,(2*j+1)**2)/W[2*j+1] for j in range(k+1))
    print(n, p[n], W[n]*s==p[n], " c_j/(2j+1) form:", W[n]*sum(F(comb(2*j,j),4**j*(2*j+1)) for j in range(k+1))==p[n])
# starting value for kernel 1/(1-ts): int theta^2 log((1+sin)/(1-sin)) dtheta
S=m.quad(lambda t: t**2*(m.log(1+m.sin(t))-m.log(2*m.sin((pi/2-t)/2)**2)),[0,pi/4,pi/2])
names=['beta4','ImLi4','ImLi3*log2','G*log2^2','pi^2*G','pi^3*log2','pi*log2^3','pi*zeta3']
vals=[b4,I4,I3*L2,G*L2**2,pi**2*G,pi**3*L2,pi*L2**3,pi*z3]
rel=m.pslq([S]+vals,maxcoeff=10**6,maxsteps=10**6)
print("S'(0) =", m.nstr(S,20), " relation:", None if rel is None else " + ".join(f"{m.mpf(-c)/rel[0]}*{n}" for c,n in zip(rel[1:],names) if c))
