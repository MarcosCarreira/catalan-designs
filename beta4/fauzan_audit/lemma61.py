import mpmath as mp
mp.mp.dps=40
T=[(3906748086,8992695531,10515596180),(2312248264,15340997855,29471737793),(1402286665,25730180724,42934365099),
(881725356,41909578246,58204231966),(578197906,65851089563,69037621310),(396324613,99481037884,78873099189),
(283911191,144325727458,84856120711),(212206188,201105762729,88396082127),(165097686,269345996903,88303382125),
(133347132,347089554156,85472321255),(111522114,430806704415,78899184238),(96349355,515561896511,70353471918),
(85815639,595448778546,58838976615),(78667711,664241383483,44421321106),(74129565,716160577112,30462865791),
(71741310,746637295669,5959622577)]
S=mp.mpf(10)**12
A=[mp.mpf(a)/S for a,b,c in T]; B=[mp.mpf(b)/S for a,b,c in T]; C=[mp.mpf(c)/S for a,b,c in T]
from fractions import Fraction as F
print('sum c =',sum(F(c,10**12) for a,b,c in T),' nested:',all(A[i]>A[i+1] and B[i]<B[i+1] for i in range(15)), ' min len', min(B[i]-A[i] for i in range(16))>mp.mpf(1)/225, 'b16<2',B[15]<2)
al=mp.mpf(3)/40; lam=mp.mpf(37)/40
def Uarc(t,a,b):
    if a<=t<=b: return mp.log((b-a)/4)
    m=(a+b)/2; return mp.log((abs(t-m)+mp.sqrt((t-a)*(t-b)))/2)
def U(t): return sum(C[j]*Uarc(t,A[j],B[j]) for j in range(16))
def V(t):
    st=mp.sqrt(t)
    return mp.log(1+t)-6*al*mp.log(t+al**2)-2+12*al+2*st*(mp.pi+mp.atan(1/st)-6*mp.atan(al/st)) if t>0 else -12*al*mp.log(al)-2+12*al
M0=mp.mpf(-1329)/200
# dense scan on [0, 4] with geometric refinement near 0
pts=[mp.mpf(10)**(-8+8*k/2000) for k in range(2001)]+[mp.mpf(k)/20000 for k in range(1,80001)]
worst=max((2*U(t)-V(t),t) for t in pts)
print('max 2U-V on grid (0,4]:',mp.nstr(worst[0],12),'at t =',mp.nstr(worst[1],8),' M0 =',float(M0),' ok:',worst[0]<=M0)
# I(rho) by (A.2) and the constant
Sj=0; I=0
for j in range(16):
    S0=Sj; Sj+=C[j]; I+=(Sj**2-S0**2)*mp.log((B[j]-A[j])/4)
Cs=-2*lam+12*al*lam*(1-mp.log(al))+3*lam**2-2*lam**2*mp.log(2*lam)
val=lam*M0-I+Cs
print('I(rho) =',mp.nstr(I,15),' C* =',mp.nstr(Cs,15),' lam*M0 - I + C* =',mp.nstr(val,15),' U =',-2733991/2000000,' ok:',val<=mp.mpf(-2733991)/2000000)
