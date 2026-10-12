"""Exact arrays for the measure nu = |t| (arcsin^2 t - pi |arcsin t|) dt / sqrt(1-t^2) and for mu = |t| dt / sqrt(1-t^2).
Every value is a dict {'1': rational, 'b4': ..., 'G': ...}.  Indices: M(i,j) = int dnu(t) t^i int_0^1 s^j/(1-ts) ds.
Facts used (route 3 notes): nu_{2k} = -2 m_{2k} O2(k), O2(k) = sum_{l<=k} (2l+1)^-2;
K-_d (d >= 2): d K_d = (d-1) K_{d-2} - sum_{v<d-2} nu_v/(v+3) - 2 (KA_d - 4G + P1_d)/d - 2 P2_d,
  KA_d = km[d] + 4G c(d/2) [d even], P1_d = sum_{n<d} m_n/(n+1), P2_d = sum_{n<d} m_n/(n+1)^2;
K-_0 = -8 beta4, K-_1 = 8G - 12 (identified by PSLQ at 60-200 digits).
M(i,j) for i >= j: K-_{i-j} - sum_{l=1}^{j} nu_{i-j+l-1}/l."""
from fractions import Fraction as Fr
from math import comb
def c(l): return Fr(comb(2*l, l), 4**l)
def add(a, b, f=1):
    out = dict(a)
    for k, v in b.items(): out[k] = out.get(k, Fr(0)) + f*v
    return out
def build_nu(top):
    m = [Fr(0)]*(top+3); m[0] = Fr(2)
    for i in range(2, top+3, 2): m[i] = i*m[i-2]/(i+1)
    km = [Fr(0)]*(top+3); km[1] = Fr(2)
    for d in range(2, top+3): km[d] = ((d-1)*km[d-2] + m[d-2] + m[d-1])/d
    nu = [Fr(0)]*(top+3); O2 = Fr(0)
    for k in range(0, (top+2)//2+1):
        O2 += Fr(1, (2*k+1)**2)
        if 2*k < len(nu): nu[2*k] = -2*m[2*k]*O2
    K = [None]*(top+1)
    K[0] = {'b4': Fr(-8)}; K[1] = {'1': Fr(-12), 'G': Fr(8)}
    P1 = [Fr(0)]*(top+2); P2 = [Fr(0)]*(top+2)
    for d in range(1, top+2):
        P1[d] = P1[d-1] + m[d-1]/d; P2[d] = P2[d-1] + m[d-1]/d**2
    S3 = [Fr(0)]*(top+2)   # S3[d] = sum_{v<d-2} nu_v/(v+3)
    for d in range(3, top+2): S3[d] = S3[d-1] + nu[d-3]/d
    for d in range(2, top+1):
        KA = {'1': km[d]}
        if d % 2 == 0: KA['G'] = 4*c(d//2)
        rhs = {k: (d-1)*v for k, v in K[d-2].items()}
        rhs = add(rhs, {'1': -S3[d] - 2*P2[d]})
        inner = add(KA, {'G': Fr(-4), '1': P1[d]})
        rhs = add(rhs, inner, Fr(-2, d))
        K[d] = {k: v/d for k, v in rhs.items() if v != 0}
    return m, km, nu, K
def M_nu(i, j, nu, K):
    assert i >= j
    out = dict(K[i-j]); out['1'] = out.get('1', Fr(0)) - sum((nu[i-j+l-1]/l for l in range(1, j+1)), Fr(0))
    return out
