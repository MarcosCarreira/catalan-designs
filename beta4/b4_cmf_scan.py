"""b4_cmf_scan.py -- recurrence (CMF) scan for beta(4), with calibrations.

Run:  python3 b4_cmf_scan.py        (no arguments; every switch is set below)
Needs only mpmath and the standard library; cmf_lite.py, rz_h1.py and rz_scout.py in the same folder.
(cmf_lite.py is a dependency-free implementation of the pFq field of ramanujantools 0.0.1; gate C0 checks it against
values computed with ramanujantools in the sandbox on 8 October.)
Writes b4_cmf_scan_log_<stamp>.txt and b4_cmf_scan_out_<stamp>.json next to this file.
Expected runtime with 2 workers: about 10-20 minutes (section A dominates).

Exponent (as ramanujantools CMF.delta): |p_n/q_n - L| = q_n^-(1+delta) with p_n/q_n reduced, L from depth 2n;
delta > 0 along a trajectory proves L irrational (Apery: delta about 0.08 for zeta(3)).
Trajectories that meet a singular unit matrix on every one of 202 monotone paths in some step are skipped and counted.

C  calibrations: C0 cmf_lite against ramanujantools, 4 trajectories at depth 40, |difference| < 1e-7 (gate);
   C1 zeta(3), Apery's recurrence, delta(100) in [0.05, 0.20] (gate; 0.135 in the sandbox);
   C2 Catalan, 3F2(1/2,1/2,1; 3/2,3/2; -1), all trajectories with steps 0..2 (x0, x1 as a multiset), delta at depths 50, 100.
A  beta(4) = 5F4(1/2,1/2,1/2,1/2,1; 3/2,3/2,3/2,3/2; -1): all trajectories with steps 0..2 (x0..x3 and y0..y3 as multisets),
   delta at depths 50 and 100; then the 30 best trajectories again from 2 shifted starts (x0..x3 + 1; y0..y3 + 1).
D  identification of the limits of the 10 best trajectories of A: p_n/q_n at depth 400 at high precision, PSLQ against
   {1, beta(4), G, pi, pi^2, pi^3, pi^4, pi log2, pi log^2 2, log2, zeta(3), Im Li3((1+i)/2), Im Li4((1+i)/2)}, with the
   number of trustworthy digits estimated from delta; a relation is reported only if PSLQ returns the same one at that
   precision and at 3/4 of it.
B  Rivoal-Zudilin pairs: very-well-poised forms S_n(eta) = c0 + c2 beta(2) + c4 beta(4) (rz_h1.py, exact), eta = (e0; e1..e5; e6)
   with 3 <= e0 <= 8, 1 <= e_j, 2 e_j < e0, sum e_j <= 5 e0 / 2 (e1..e5 as a multiset); for every pair, the beta(2)-free
   combination L_n = c2(eta') S_n(eta) - c2(eta) S_n(eta'); balance B_n = log|L_n| + log den(L_n) at n = 8 and 12 and the
   slope (B_12 - B_8)/4 (a proof needs a negative slope). Gate B0: exact form = direct series for eta = (3;1,1,1,1,1;1), n = 2.
"""
import os, sys, math, time, json, datetime, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE); sys.path.insert(0, HERE)
from fractions import Fraction as F
from multiprocessing import Pool
import mpmath as mp

# ---------------- switches ----------------
WORKERS = 2
STEPS = (0, 1, 2)
DEPTHS = (50, 100)
TOP_SHIFT = 30
TOP_ID = 10
ID_DEPTH = 400
RZ_E0 = range(3, 9)
RZ_NS = (8, 12)
TOP_PAIRS = 15
# ------------------------------------------

STAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
LOGF = os.path.join(HERE, f"b4_cmf_scan_log_{STAMP}.txt"); OUTF = os.path.join(HERE, f"b4_cmf_scan_out_{STAMP}.json")
T0 = time.time(); _log = None
def say(s):
    global _log
    line = f"[{time.time()-T0:7.1f}s] {s}"; print(line, flush=True)
    if _log is None: _log = open(LOGF, "w")
    _log.write(line + "\n"); _log.flush()
FAILS = {}
def gate(name, ok, info=""):
    FAILS[name] = FAILS.get(name, 0) + (0 if ok else 1)
    if not ok: say(f"FAIL {name} {info}")

import cmf_lite as CL
H = F(1, 2)
STARTS = {'b4': ((H, H, H, H, F(1)), (3*H, 3*H, 3*H, 3*H)), 'cat': ((H, H, F(1)), (3*H, 3*H))}
def start_of(kind, shift):
    xs, ys = STARTS[kind]; xs = list(xs); ys = list(ys)
    if shift == 'x+1': xs = [v+1 for v in xs[:-1]] + [xs[-1]]
    elif shift == 'y+1': ys = [v+1 for v in ys]
    return xs, ys
def traj_task(t):
    kind, xs, ys, shift = t
    x0, y0 = start_of(kind, shift)
    out = dict(kind=kind, xs=list(xs), ys=list(ys), shift=shift)
    try:
        out['delta'] = CL.deltas(x0, y0, -1, xs, ys, DEPTHS)
    except ArithmeticError as e:
        out['error'] = 'singular'
    except Exception as e:
        out['error'] = repr(e)[:120]
    return out
def trajectories(kind):
    nx, ny = (5, 4) if kind == 'b4' else (3, 2)
    sym_x = nx - 1                       # first nx-1 numerator parameters are equal; last one separate
    xs_list = [tuple(ms) + (last,) for ms in itertools.combinations_with_replacement(STEPS, sym_x) for last in STEPS]
    ys_list = list(itertools.combinations_with_replacement(STEPS, ny))
    return [(kind, xs, ys, 'none') for xs in xs_list for ys in ys_list if sum(xs)+sum(ys) > 0]
def good(r): return 'delta' in r and all(math.isfinite(d) for d in r['delta'])

def rz_task(eta):
    src = open('rz_h1.py').read().split('ETA1=')[0]
    g = {'__name__': 'rzmod'}; exec(src, g)
    out = {}
    for n in RZ_NS:
        (c0, c2, c4), _ = g['full_form'](eta, n); out[n] = (c0, c2, c4)
    return eta, out

if __name__ == "__main__":
    say(f"b4_cmf_scan.py  workers={WORKERS} steps={STEPS} depths={DEPTHS}")
    REF = [('b4', (1, 1, 0, 2, 1), (0, 0, 1, 2), -0.3890847581732966), ('b4', (0, 0, 0, 0, 1), (1, 1, 1, 1), -0.69421964),
           ('b4', (1, 1, 1, 1, 2), (2, 2, 2, 2), -0.80406388), ('cat', (1, 1, 1), (1, 1), -0.94998351)]
    worst = 0.0
    for kind, xs, ys, ref in REF:
        x0, y0 = start_of(kind, 'none'); d = CL.deltas(x0, y0, -1, xs, ys, [40])[0]; worst = max(worst, abs(d-ref))
    gate("C0", worst < 1e-7, f"worst {worst}")
    say(f"C0 cmf_lite vs ramanujantools reference values (depth 40): worst difference {worst:.1e}")
    d = CL.apery_delta(100)
    gate("C1", 0.05 <= d <= 0.20, f"delta {d}")
    say(f"C1 zeta(3), Apery's recurrence, delta(100) = {d:.4f} (asymptotic about 0.08)")
    # ---------------- C2 and A
    res = {}
    with Pool(WORKERS) as pool:
        for kind in ('cat', 'b4'):
            tasks = trajectories(kind); out = []
            say(f"{'C2 Catalan 3F2' if kind == 'cat' else 'A beta(4) 5F4'}: {len(tasks)} trajectories")
            for i, r in enumerate(pool.imap_unordered(traj_task, tasks, chunksize=4)):
                out.append(r)
                if (i+1) % 100 == 0: say(f"   {i+1}/{len(tasks)} done")
            ok = [r for r in out if good(r)]
            ok.sort(key=lambda r: -r['delta'][-1])
            errs = len(out) - len(ok); sing = sum(1 for r in out if r.get('error') == 'singular')
            ds = sorted(r['delta'][-1] for r in ok)
            say(f"   {len(ok)} valid, {errs} not valid ({sing} skipped as singular); median delta({DEPTHS[-1]}) = {ds[len(ds)//2]:.3f}; best ten:")
            for r in ok[:10]: say(f"     x {r['xs']} y {r['ys']}: delta {['%.4f' % x for x in r['delta']]}")
            res[kind] = out
        best = [r for r in res['b4'] if good(r)]; best.sort(key=lambda r: -r['delta'][-1]); best = best[:TOP_SHIFT]
        tasks = [('b4', tuple(r['xs']), tuple(r['ys']), s) for r in best for s in ('x+1', 'y+1')]
        say(f"A2 shifted starts for the {len(best)} best beta(4) trajectories")
        sh = list(pool.imap_unordered(traj_task, tasks))
        shok = [r for r in sh if good(r)]; shok.sort(key=lambda r: -r['delta'][-1])
        for r in shok[:10]: say(f"     {r['shift']} x {r['xs']} y {r['ys']}: delta {['%.4f' % x for x in r['delta']]}")
        res['b4_shift'] = sh
    # ---------------- D identification
    mp.mp.dps = 400
    pi = mp.pi; L2 = mp.log(2); w = mp.mpc(1, 1)/2
    basis = {'1': 1, 'beta4': (mp.zeta(4, mp.mpf(1)/4)-mp.zeta(4, mp.mpf(3)/4))/256, 'G': +mp.catalan, 'pi': pi, 'pi^2': pi**2,
             'pi^3': pi**3, 'pi^4': pi**4, 'pi*log2': pi*L2, 'pi*log2^2': pi*L2**2, 'log2': L2, 'zeta3': mp.zeta(3),
             'ImLi3': mp.im(mp.polylog(3, w)), 'ImLi4': mp.im(mp.polylog(4, w))}
    allb = [r for r in res['b4'] + res['b4_shift'] if good(r)]; allb.sort(key=lambda r: -r['delta'][-1])
    ident = []
    for r in allb[:TOP_ID]:
        x0, y0 = start_of('b4', r['shift'])
        try:
            p, q = CL.walk(x0, y0, -1, r['xs'], r['ys'], [ID_DEPTH])[ID_DEPTH]
            g_ = math.gcd(p, q); p //= g_; q //= g_
            digits = int(0.9*(1+r['delta'][-1])*math.log10(abs(q))) if q not in (0, 1, -1) else 0
            mp.mp.dps = 400
            L = mp.mpf(p)/q
            names = list(basis); rels = []
            for dd in (min(digits, 380), int(0.75*min(digits, 380))):   # a relation counts only if found at both precisions
                mp.mp.dps = max(30, dd)
                rels.append(mp.pslq([L]+[+basis[k] for k in names], maxcoeff=10**8, maxsteps=10**6) if dd >= 60 else None)
            mp.mp.dps = 400
            rel = rels[0] if (rels[0] is not None and rels[0] == rels[1]) else None
            idn = None if rel is None else {k: str(F(-c, rel[0])) for c, k in zip(rel[1:], names) if c}
            ident.append(dict(traj=r, digits=digits, L=mp.nstr(L, 25), relation=idn))
            say(f"D  x {r['xs']} y {r['ys']} {r['shift']}: L = {mp.nstr(L, 20)}, about {digits} digits; relation {idn}")
        except Exception as e:
            say(f"D  x {r['xs']} y {r['ys']} {r['shift']}: error {repr(e)[:100]}")
    # ---------------- B Rivoal-Zudilin pairs
    src = open('rz_h1.py').read().split('ETA1=')[0]
    g = {'__name__': 'rzmod'}; exec(src, g)
    mp.mp.dps = 60
    (c0, c2, c4), _ = g['full_form']((3, 1, 1, 1, 1, 1, 1), 2)
    b2v = +mp.catalan; b4v = (mp.zeta(4, mp.mpf(1)/4)-mp.zeta(4, mp.mpf(3)/4))/256
    fv = lambda x: mp.mpf(x.numerator)/x.denominator
    exact = fv(c0)+fv(c2)*b2v+fv(c4)*b4v; direct = g['direct']((3, 1, 1, 1, 1, 1, 1), 2)
    gate("B0", abs(exact-direct) < mp.mpf('1e-40'), f"diff {mp.nstr(exact-direct, 3)}")
    say(f"B0 Rivoal-Zudilin exact form vs direct series: diff {mp.nstr(exact-direct, 3)}")
    etas = []
    for e0 in RZ_E0:
        vals = [e for e in range(1, e0) if 2*e < e0]
        for ms in itertools.combinations_with_replacement(vals, 5):
            for e6 in vals:
                if 2*(sum(ms)+e6) <= 5*e0: etas.append((e0,)+ms+(e6,))
    say(f"B  {len(etas)} parameter vectors eta; forms at n = {RZ_NS}")
    forms = {}
    with Pool(WORKERS) as pool:
        for eta, out in pool.imap_unordered(rz_task, etas):
            forms[eta] = out
    mp.mp.dps = 120
    b4v = (mp.zeta(4, mp.mpf(1)/4)-mp.zeta(4, mp.mpf(3)/4))/256
    pairs = []
    keys = sorted(forms)
    for i in range(len(keys)):
        for j in range(i+1, len(keys)):
            B = {}
            for n in RZ_NS:
                a0, a2, a4 = forms[keys[i]][n]; e0_, e2, e4 = forms[keys[j]][n]
                A = e2*a4 - a2*e4; Bc = e2*a0 - a2*e0_
                if A == 0: B = None; break
                val = fv(A)*b4v + fv(Bc)
                den = A.denominator*Bc.denominator//math.gcd(A.denominator, Bc.denominator)
                B[n] = float(mp.log(abs(val))) + math.log(den) if val != 0 else float('-inf')
            if B: pairs.append(dict(eta1=keys[i], eta2=keys[j], B=B, slope=(B[RZ_NS[1]]-B[RZ_NS[0]])/(RZ_NS[1]-RZ_NS[0])))
    pairs.sort(key=lambda r: r['slope'])
    sl = sorted(p['slope'] for p in pairs)
    say(f"B  {len(pairs)} pairs; slope per n: min {sl[0]:.3f}, median {sl[len(sl)//2]:.3f}; best {TOP_PAIRS}:")
    for p in pairs[:TOP_PAIRS]: say(f"     {p['eta1']} & {p['eta2']}: B_{RZ_NS[0]} {p['B'][RZ_NS[0]]:+.1f}, B_{RZ_NS[1]} {p['B'][RZ_NS[1]]:+.1f}, slope {p['slope']:+.3f}")
    total = sum(FAILS.values())
    json.dump(dict(stamp=STAMP, gates=FAILS, C2=res['cat'], A=res['b4'], A2=res['b4_shift'], D=ident,
                   B=pairs[:200], gate_fails=total), open(OUTF, "w"), indent=1, default=str)
    say(f"gates: {FAILS}")
    say(f"GATE_FAILS = {total}")
    say(f"log saved to {LOGF}")
    say(f"json saved to {OUTF}")
