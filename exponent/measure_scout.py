"""measure_scout.py (scout) -- how far from G the nearest zero of the Catalan determinant lies.

For the preprint's determinant of size n (design ratio 48:7:4:4:2, rounded as in full_cat.py), the matrix is
Frat + X*FG with X in place of G, so Delta(X) is a polynomial of degree <= n. Fauzan's Appendix C route to an
irrationality measure needs |Delta(x)| <= 2|Delta(G)| for all x within e^(-rho*n) of G, with rho bounded.
The scout expands P(G+d) = sum_k c_k d^k for the primitive integer polynomial P of Delta and prints the
largest d* with sum_k |c_k| d*^k <= 2|P(G)|, as -log d*/n (bounded for a power measure) and -log d*/n^2,
and Lambda = log|P(G)|/n^2.
Usage: python3 measure_scout.py [n1,n2,...]   (default SIZES with RATIO, both set at the top; writes measure_scout_log_<time>.txt)
Also prints the Hadamard constant: (1/n^2) log of the product of row lengths of [A | FG], a bound for all mixed minors.
Requirements: python-flint, mpmath. Sandbox timing: n=24 1 s, n=48 11 s; n=72 takes over a minute, larger n longer."""
import sys, math, time
from fractions import Fraction as Fr
import mpmath as mpm
import flint
sys.path.insert(0, '.')
from full_cat import build
try: sys.set_int_max_str_digits(0)
except AttributeError: pass
T0 = time.time()
# Design ratio (n : b : q : g : h). OURS = Theorem design (96,15,6,8,5); OPENAI = preprint (48,7,4,4,2).
OURS, OPENAI = (96, 15, 6, 8, 5), (48, 7, 4, 4, 2)
RATIO = OURS
SIZES = [192]      # 48-144 done in log _192255; exact ratio at 96 and 192, other sizes rounded
STAMP = time.strftime('%Y%m%d_%H%M%S')
LOG = f'measure_scout_log_{STAMP}.txt'
def say(msg):
    print(msg, flush=True)
    with open(LOG, 'a') as fh: fh.write(msg + '\n')

def delta_poly(n, b, q, g, h):
    Frat, FG = build(n, b, q, g, h)
    A0, A1 = [], []
    for r in range(n):
        L = 1
        for x in Frat[r] + FG[r]: L = L*x.denominator//math.gcd(L, x.denominator)
        A0.append([int(x*L) for x in Frat[r]]); A1.append([int(x*L) for x in FG[r]])
    xs = list(range(n+1)); ys = []
    for x in xs:
        M = flint.fmpz_mat([[A0[r][k] + x*A1[r][k] for k in range(n)] for r in range(n)])
        ys.append(M.det())
    # Newton divided differences, exact
    coef = [Fr(int(y)) for y in ys]
    for j in range(1, n+1):
        for i in range(n, j-1, -1):
            coef[i] = (coef[i]-coef[i-1])/(xs[i]-xs[i-j])
    P = flint.fmpq_poly([0])
    basis = flint.fmpq_poly([1])
    for j in range(n+1):
        P += flint.fmpq(coef[j].numerator, coef[j].denominator)*basis
        basis *= flint.fmpq_poly([-xs[j], 1])
    num = P.numer()
    cont = num.content() if hasattr(num, 'content') else None
    Z = flint.fmpz_poly([int(c) for c in num.coeffs()])
    gcdc = 0
    for c in Z.coeffs(): gcdc = math.gcd(gcdc, int(c))
    Z = flint.fmpz_poly([int(c)//gcdc for c in Z.coeffs()])
    return Z, Frat, FG

def taylor_at_G(Z, dps):
    """Coefficients c_k of P(G+d) = sum_k c_k d^k, by repeated synthetic division at G."""
    mpm.mp.dps = dps
    G = +mpm.catalan
    a = [mpm.mpf(int(c)) for c in reversed(Z.coeffs())]   # highest degree first
    out = []
    while a:
        acc = mpm.mpf(0); b = []
        for x in a:
            acc = acc*G + x; b.append(acc)
        out.append(b[-1]); a = b[:-1]
    return out

def radius(c):
    """Largest d with sum_k |c_k| d^k <= 2|c_0|, found by bisection in log d."""
    lc = [mpm.log(abs(x)) if x != 0 else None for x in c]
    def f(ld):
        terms = [l + k*ld for k, l in enumerate(lc) if l is not None]
        m = max(terms)
        return m + mpm.log(sum(mpm.exp(t-m) for t in terms)) - (mpm.log(2) + lc[0])
    lo, hi = mpm.mpf(-10**6), mpm.mpf(0)
    for _ in range(200):
        mid = (lo+hi)/2
        if f(mid) <= 0: lo = mid
        else: hi = mid
    return lo

def hadamard(Frat, FG, n):
    """(1/n^2) * log of prod_r (sqrt(|row of A|^2 + |row of FG|^2)), A = Frat + G*FG: a Hadamard bound
    for every minor mixing columns of A and FG, i.e. for every Taylor coefficient of Delta at G (up to binomials)."""
    G = float(mpm.catalan)
    tot = 0.0
    for r in range(n):
        sa = sum((float(Frat[r][k]) + G*float(FG[r][k]))**2 for k in range(n))
        sb = sum(float(FG[r][k])**2 for k in range(n))
        tot += 0.5*math.log(sa + sb)
    return tot/n**2

def run(n):
    R = RATIO
    d = (n,) + tuple(round(x*n/R[0]) for x in R[1:])
    t = time.time()
    Z, Frat, FG = delta_poly(*d)
    Hc = hadamard(Frat, FG, n)
    deg = Z.degree()
    digits = len(str(max(abs(int(c)) for c in Z.coeffs())))
    res = []
    base = digits + n*n//4 + 60          # |P(G)| is about 10^(-n^2/4) or larger at these sizes
    for dps in (base, base + 200):
        c = taylor_at_G(Z, dps)
        res.append((c, radius(c)))
    (c0, r0), (c1, r1) = res
    mpm.mp.dps = base + 200
    Lam = float(mpm.log(abs(c1[0])))/n**2
    rel = float(abs(r0-r1))
    ld = float(r1)
    say(f"n={n:4d} design={d}: deg P={deg}, log10 max coeff={digits}, Lambda=log|P(G)|/n^2={Lam:+.5f}; "
          f"radius d* with |P| <= 2|P(G)| on |x-G| <= d*: log d* = {ld:.2f}, -log d*/n = {-ld/n:.4f}, -log d*/n^2 = {-ld/n**2:.5f} "
          f"(prec diff {rel:.1e}); Hadamard constant {Hc:.4f}  [{time.time()-t:.1f}s, total {time.time()-T0:.0f}s]")

if __name__ == "__main__":
    ns = list(map(int, sys.argv[1].split(','))) if len(sys.argv) > 1 else SIZES
    say(f'measure_scout.py  ratio {RATIO}  sizes {ns}')
    for n in ns: run(n)
    say(f'log saved to {LOG}')
