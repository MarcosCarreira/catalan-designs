import mpmath as mp, math
mp.mp.dps=30
def G(y):
    y=mp.mpf(y)
    if y==0: return 5*mp.log(4)-8*mp.log(64)+12*mp.log(6)-6*mp.log(2)
    return 5*mp.log(y**2+4)-8*mp.log(y**2+64)+5*y*mp.atan(2/y)-2*y*mp.atan(8/y)+12*mp.log(6)-6*mp.log(2)
lam0=mp.diff(G,mp.mpf(4)/5); print('lambda0 = G\'(4/5) =',mp.nstr(lam0,10),' paper 3.009194; formula 5atan(5/2)-2atan(10) =',mp.nstr(5*mp.atan(2.5)-2*mp.atan(10),10))
c1=G(mp.mpf(4)/5)-mp.mpf(4)/5*lam0; c2=2*mp.pi-lam0
c3=9*mp.log(2*c2/mp.pi)-3*c1
print('c1 =',mp.nstr(c1,10),'(<0)',' c2 =',mp.nstr(c2,10),'(>0)',' c3 =',mp.nstr(c3,10),' paper 31.614437')
c3p=9*mp.log(1+(10*mp.atan(mp.mpf(2)/5)-4*mp.atan(mp.mpf(1)/10))/mp.pi)-3*(12*mp.log(3)-32*mp.log(2)+5*mp.log(mp.mpf(29)/25)-8*mp.log(mp.mpf(101)/100))
print('c3 by the paper closed form =',mp.nstr(c3p,10))
print('max over y in [0,200] of G(y)-lambda0*y - c1:',mp.nstr(max(G(mp.mpf(k)/100)-lam0*mp.mpf(k)/100-c1 for k in range(0,20001)),6),'(must be <= 0)')
# Lemma 6.2: |R_n(2n+1/2+iny)| <= (10/n) e^{nG(y)}
def logR(n,z):
    s=2*mp.log(mp.factorial(6*n))-3*mp.log(mp.factorial(2*n))
    s+=3*sum(mp.log(abs(z-2*n+k)) for k in range(2*n))-2*sum(mp.log(abs(z+k)) for k in range(6*n+1))
    return s
worst=-1e9
for n in (1,2,5,10,20,50):
    for k in range(0,401):
        y=mp.mpf(k)/20; z=2*n+mp.mpf(1)/2+1j*n*y
        d=logR(n,z)-(mp.log(mp.mpf(10)/n)+n*G(y)); worst=max(worst,d)
print('Lemma 6.2: max of log|R_n| - log((10/n)e^{nG}) over n in {1..50}, y in [0,20]:',mp.nstr(worst,6),'(must be <= 0)')
# Lemma 6.5 core identity (6.9) for N = 1..4 via Hankel determinants of the moments of 2/(e^{pi y}+e^{-pi y})
for N in range(1,5):
    mom=[mp.quad(lambda y:y**k/mp.cosh(mp.pi*y),[-mp.inf,0,mp.inf]) for k in range(2*N-1)]
    H=mp.det(mp.matrix([[mom[i+j] for j in range(N)] for i in range(N)]))
    print(f'(6.9) N={N}: Hankel det of moments = {mp.nstr(H,15)}  prod k!^2/4^k = {mp.nstr(mp.fprod(mp.factorial(k)**2/4**k for k in range(N)),15)}')
# Lemma 6.1 at n = 1: contour integral vs series
n=1; a,b=2,-3
from fractions import Fraction as Fr
def R(z): return mp.factorial(6*n)**2/mp.factorial(2*n)**3*mp.rf(z-2*n,2*n)**3/mp.rf(z,6*n+1)**2
for s in (0,1,2):
    f=lambda x: R(2*n+mp.mpf(1)/2+1j*x)*(2*n+mp.mpf(1)/2+1j*x)**s*(mp.pi/mp.sin(mp.pi*(2*n+mp.mpf(1)/2+1j*x)))**2*(3*a-2*b*mp.pi*mp.cot(mp.pi*(2*n+mp.mpf(1)/2+1j*x)))
    I=mp.quad(f,[-mp.inf,-5,0,5,mp.inf])/(2*mp.pi)
    g=lambda t: R(t)*t**s
    ser=-3*a*mp.nsum(lambda m: mp.diff(g,m),[1,mp.inf])+b*mp.nsum(lambda m: mp.diff(g,m,2),[1,mp.inf])
    print(f'Lemma 6.1 n=1, P=t^{s}: contour {mp.nstr(I,15)}  series {mp.nstr(ser,15)}')
