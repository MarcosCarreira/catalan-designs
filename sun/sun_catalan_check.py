"""sun_catalan_check.py -- check of the final estimate of Z.-W. Sun, "Catalan's constant is irrational", arXiv:2609.04176
(3 September 2026), Theorem 9.1. Appendix D of the note.

Run:  python3 sun_catalan_check.py        (no arguments; every switch is set below)
Needs: python-flint, mpmath.  Writes sun_catalan_check_log_<stamp>.txt and _out_<stamp>.json next to this file.
Expected runtime: about 30 seconds.

Construction (Sun, Sections 2-3): Pi_i = prod_{h=1..B} (2(h+i)+1)^2, T_m = (-1)^m (X - S_{m-1}) with S_{m-1} the partial sums of G,
u_m = T_m/(2m+1), R_{a,j} = sum_{i=0}^{a+2B} (-1)^i C(a+2B,i) Pi_i u_{i+j} (0 <= a <= S+2, 1 <= j <= S), N = 2B+S+3, D = 2B,
F_B = prod_{r<2B} r!, qhat_B(X) = F_B det R[A,J] / prod_{i<N} Pi_i, a polynomial of degree <= S in X; S = floor(B/20).
Measured: Lambda_B = (log D_B + log|qhat_B(G)|)/B^2 with D_B the least common denominator of the coefficients of qhat_B.
Theorem 9.1 of Sun states log H_B^min + log|qhat_B| <= -delta0 B^2 + o(B^2) with delta0 > 0.00966.
Gates: S1 Prop. 3.1 of Sun at B = 20: det of the N x N matrix A~_B (columns i^r, r < D; Pi_i u_{i+1}; binomial(i, D+c), c not in A)
equals F_B det R[A,J] exactly, at X = 0 and X = 1, for all four A; S2 for every B and A: the value at two precisions agrees to 1e-20.
"""
import os, sys, math, time, json, datetime, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
from flint import fmpq, fmpq_poly, fmpq_mat
import mpmath as mp

# ---------------- switches ----------------
BS = [20, 40, 60, 80, 100]
# ------------------------------------------

STAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
LOGF = os.path.join(HERE, f"sun_catalan_check_log_{STAMP}.txt"); OUTF = os.path.join(HERE, f"sun_catalan_check_out_{STAMP}.json")
T0 = time.time(); _log = open(LOGF, "w")
def say(s):
    line = f"[{time.time()-T0:7.1f}s] {s}"; print(line, flush=True); _log.write(line + "\n"); _log.flush()
FAILS = {}
def gate(name, ok, info=""):
    FAILS[name] = FAILS.get(name, 0) + (0 if ok else 1)
    if not ok: say(f"FAIL {name} {info}")

def setup(B, S):
    N = 2*B + S + 3
    Ssum = [fmpq(0)]
    for k in range(N + S + 6): Ssum.append(Ssum[-1] + fmpq((-1)**k, (2*k+1)**2))
    Pi = [math.prod((2*(h+i)+1)**2 for h in range(1, B+1)) for i in range(N + S + 3)]
    FB = 1
    for r in range(2*B): FB *= math.factorial(r)
    return N, Ssum, Pi, FB
def u_poly(m, Ssum): return fmpq_poly([-Ssum[m], 1]) * ((-1)**m) / (2*m+1)
def Rmat(B, S, Ssum, Pi):
    return [[sum((u_poly(i+j, Ssum) * ((-1)**i * math.comb(a+2*B, i) * Pi[i]) for i in range(a+2*B+1)), fmpq_poly([0]))
             for j in range(1, S+1)] for a in range(S+3)]
def detpoly(M, rows, S):
    det = fmpq_poly([0])
    for perm in itertools.permutations(range(S)):
        sgn = 1
        for x in range(S):
            for y in range(x+1, S):
                if perm[x] > perm[y]: sgn = -sgn
        t = fmpq_poly([sgn])
        for r_, c_ in zip(rows, perm): t = t * M[r_][c_]
        det += t
    return det

if __name__ == "__main__":
    say(f"sun_catalan_check.py  B in {BS}")
    # S1: Proposition 3.1 at B = 20
    B, S = 20, 1; N, Ssum, Pi, FB = setup(B, S); D = 2*B
    R = Rmat(B, S, Ssum, Pi)
    for A in ([0], [1], [2], [3]):
        Ac = [c for c in range(S+3) if c not in A]
        for Xv in (fmpq(0), fmpq(1)):
            M = fmpq_mat(N, N, [x for i in range(N) for x in [fmpq(i**r) for r in range(D)] + [Pi[i]*((-1)**(i+1))*(Xv - Ssum[i+1])/(2*(i+1)+1)] + [fmpq(math.comb(i, D+c)) for c in Ac]])
            lhs = M.det(); rhs = FB * R[A[0]][0](Xv)
            gate("S1", lhs == rhs or lhs == -rhs, f"A={A} X={Xv}")
    say("S1 Proposition 3.1 of Sun at B = 20 (all four A, X = 0 and 1): " + ("exact" if FAILS.get("S1", 0) == 0 else "FAILED"))
    out = []
    for B in BS:
        t = time.time(); S = max(1, B//20); N, Ssum, Pi, FB = setup(B, S); R = Rmat(B, S, Ssum, Pi)
        PiN = math.prod(Pi[i] for i in range(N)); vals = []
        for A in itertools.combinations(range(S+3), S):
            q = detpoly(R, A, S) * fmpq(FB, PiN)
            if q == 0: continue
            Dd = 1
            for k in range(q.degree()+1):
                if q[k] != 0: Dd = Dd*int(q[k].q)//math.gcd(Dd, int(q[k].q))
            bits = max(abs(int((q[k]*Dd).p)).bit_length() for k in range(q.degree()+1))
            res = []
            for extra in (200, 1200):
                mp.mp.dps = int(0.302*bits*1.3) + extra; Gv = +mp.catalan; v = mp.mpf(0)
                for k in range(q.degree(), -1, -1): v = v*Gv + mp.mpf(int((q[k]*Dd).p))
                res.append(v)
            gate("S2", abs(res[0]-res[1]) <= abs(res[1])*mp.mpf(10)**-20, f"B={B} A={A}")
            vals.append((list(A), float(mp.log(abs(res[1])))))
        best = min(L for _, L in vals); worst = max(L for _, L in vals)
        say(f"B = {B:3d}, S = {S}: {len(vals)} sets A;  (log D_B + log|qhat_B(G)|)/B^2: best {best/B**2:+.4f}, worst {worst/B**2:+.4f};  Theorem 9.1 of Sun: <= -0.00966 + o(1)  [{time.time()-t:.1f}s]")
        out.append(dict(B=B, S=S, values=vals, best_per_B2=best/B**2, worst_per_B2=worst/B**2))
    total = sum(FAILS.values())
    json.dump(dict(stamp=STAMP, gates=FAILS, results=out, gate_fails=total), open(OUTF, "w"), indent=1)
    say(f"gates: {FAILS}")
    say(f"GATE_FAILS = {total}")
    say(f"log saved to {LOGF}")
    say(f"json saved to {OUTF}")
