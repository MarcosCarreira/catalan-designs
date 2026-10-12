import mpmath as mp
mp.mp.dps=30
exec(open('lemma61.py').read().split('M0=mp.mpf')[0])
M0=mp.mpf(-1329)/200
# (6.8): for t >= 2, 2U - V <= 13/10 log t + 6 al^3/t - 2 pi sqrt t ; and RHS + sqrt(t)/K decreasing on [2,inf) with value at 2 < M0 (K>=2)
worst=-10
for k in range(0,4000):
    t=2*mp.mpf(1.003)**k
    d=(2*U(t)-V(t))-(mp.mpf(13)/10*mp.log(t)+6*al**3/t-2*mp.pi*mp.sqrt(t))
    worst=max(worst,d)
print('(6.8): max of (2U-V) - RHS over t in [2, 2e5]:',mp.nstr(worst,6),'(must be <= 0)')
rhs=lambda t,K: mp.mpf(13)/10*mp.log(t)+6*al**3/t-2*mp.pi*mp.sqrt(t)+mp.sqrt(t)/K
print('(6.7)-(6.8): RHS+sqrt(t)/K at t=2, K=2:',mp.nstr(rhs(2,2),8),' M0 =',float(M0),' derivative sign on [2,50] (K=2) max:',mp.nstr(max(mp.diff(lambda x:rhs(x,2),2+k*0.05) for k in range(960)),5))
# (6.11): w(y) <= 8192 (1+y)^5 e^{-2 pi y}
w=lambda y:(2*mp.pi)**4*y**5/12*mp.polylog(-4,mp.e**(-2*mp.pi*y))
print('(6.11): max ratio w/(8192(1+y)^5 e^{-2pi y}) on (0,60]:',mp.nstr(max(w(mp.mpf(k)/100)/(8192*(1+mp.mpf(k)/100)**5*mp.e**(-2*mp.pi*mp.mpf(k)/100)) for k in range(1,6001)),6))
# (6.12): 0 <= sum_{j<=m} log(t+(j/K)^2) - K int_0^{m/K} log(t+u^2) du <= 2 log K + 2
def I(c,t):  # int_0^c log(t+u^2) du
    if t==0: return c*mp.log(c**2)-2*c
    return c*mp.log(t+c**2)-2*c+2*mp.sqrt(t)*mp.atan(c/mp.sqrt(t))
bad=0; mx=-1
for K in (40,80,400):
    for m in (1,K//7,K//2,K):
        for t in (0,mp.mpf('1e-6'),mp.mpf('0.01'),mp.mpf('0.3'),1,5):
            dlt=sum(mp.log(t+(mp.mpf(j)/K)**2) for j in range(1,m+1))-K*I(mp.mpf(m)/K,t)
            r=dlt/(2*mp.log(K)+2); mx=max(mx,r)
            if dlt< -1e-20 or r>1: bad+=1
print('(6.12): violations',bad,' max ratio to 2logK+2:',mp.nstr(mx,5))
# 652
print('int t^-1/2 (1+sqrt t)^5 e^-sqrt t =',mp.nstr(mp.quad(lambda t:t**-0.5*(1+mp.sqrt(t))**5*mp.e**(-mp.sqrt(t)),[0,1,10,100,mp.inf]),10))
# (6.15) and C*
lam=mp.mpf(37)/40; Cs=-2*lam+12*al*lam*(1-mp.log(al))+3*lam**2-2*lam**2*mp.log(2*lam)
for n in (1,10,100,1000,10000):
    K=40*n; N=3*n; h=37*n
    logS=2*h*mp.loggamma(K+1)+(h-1)*mp.log(4)-12*h*mp.loggamma(N+1)-2*sum(mp.loggamma(2*i+1) for i in range(1,h)) if n<=1000 else None
    if logS is None: continue
    bound=(2*lam-12*al*lam-2*lam**2)*K**2*mp.log(K)+Cs*K**2+6*h*mp.log(K)+6*h
    print(f'n={n}: log S_K = {mp.nstr(logS,10)}, bound (6.15) = {mp.nstr(bound,10)}, ok={logS<=bound}, (logS - coef K^2 logK)/K^2 = {mp.nstr((logS-(2*lam-12*al*lam-2*lam**2)*K**2*mp.log(K))/K**2,8)}  C* = {mp.nstr(Cs,8)}')
