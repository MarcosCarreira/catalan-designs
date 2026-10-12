"""b4_cmf_filter.py -- what the trajectories of b4_cmf_scan.py converge to.

Run:  python3 b4_cmf_filter.py        (no arguments; every switch is set below)
Needs: mpmath; cmf_lite.py and the scan output SCAN_JSON in the same folder.
Writes b4_cmf_filter_log_<stamp>.txt and b4_cmf_filter_out_<stamp>.json next to this file.
Expected runtime with 2 workers: about 5-15 minutes.

A walk of the pFq field converges to a ratio of two linear forms, L = (sum a_i b_i)/(sum c_i b_i), b_i in a basis of constants.
For every valid trajectory of the scan: p/q at depth DEPTH, number of trustworthy digits from the scan's delta(100),
then PSLQ on (b_i, L b_i). A relation is accepted only if PSLQ returns the same one at that precision and at 3/4 of it,
and only if at least MIN_DIGITS digits are available. First a basis of pi-type constants is tried, then the full basis.
Classes: 'beta4' (beta(4) in the relation), 'G' (Catalan's constant, no beta(4)), 'pi-type' (only 1, pi, pi^2, pi^3,
log 2, pi log 2, zeta(3)), 'unidentified'. Ranking of delta(100) inside each class.
Gates: F1 the best beta(4) trajectory x(0,1,1,1,1) y(0,0,0,2) gives L = (pi-2)/(8pi-12); F2 the best Catalan trajectory
x(1,1,1) y(0,2) gives L = (pi-2)/(4pi-4); F3 Apery's recurrence: the limit is identified as zeta(3) with the full basis.
"""
import os, sys, math, time, json, datetime
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE); sys.path.insert(0, HERE)
from fractions import Fraction as F
from multiprocessing import Pool
import mpmath as mp
import cmf_lite as CL

# ---------------- switches ----------------
WORKERS = 2
SCAN_JSON = "b4_cmf_scan_out_20261008_111915.json"
DEPTH = 200
MIN_DIGITS = 120
MAX_DPS = 450
MAXCOEFF = 10**9
LIMIT_ITEMS = None      # None = all trajectories
# ------------------------------------------

STAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
LOGF = os.path.join(HERE, f"b4_cmf_filter_log_{STAMP}.txt"); OUTF = os.path.join(HERE, f"b4_cmf_filter_out_{STAMP}.json")
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

H = F(1, 2)
STARTS = {'b4': ((H, H, H, H, F(1)), (3*H, 3*H, 3*H, 3*H)), 'cat': ((H, H, F(1)), (3*H, 3*H))}
PI_NAMES = ['1', 'pi', 'pi^2', 'pi^3', 'log2', 'pi*log2', 'zeta3']
FULL = {'b4': ['1', 'beta4', 'G', 'pi', 'pi^2', 'pi^3', 'log2', 'pi*log2', 'zeta3'],
        'cat': ['1', 'G', 'pi', 'pi^2', 'log2', 'pi*log2'],
        'z3': ['1', 'zeta3', 'pi^2', 'log2', 'G']}
def consts():
    pi = mp.pi; L2 = mp.log(2)
    return {'1': mp.mpf(1), 'beta4': (mp.zeta(4, mp.mpf(1)/4)-mp.zeta(4, mp.mpf(3)/4))/256, 'G': +mp.catalan, 'pi': pi,
            'pi^2': pi**2, 'pi^3': pi**3, 'log2': L2, 'pi*log2': pi*L2, 'zeta3': mp.zeta(3)}
def moebius(p, q, names, digits):
    """stable relation sum a_i b_i - L sum c_i b_i = 0, or None"""
    rels = []
    for dd in (min(digits, MAX_DPS), int(0.75*min(digits, MAX_DPS))):
        mp.mp.dps = dd
        C = consts(); L = mp.mpf(p)/q
        vec = [C[k] for k in names] + [L*C[k] for k in names]
        try: rels.append(mp.pslq(vec, maxcoeff=MAXCOEFF, maxsteps=10**6))
        except Exception: rels.append(None)
    if rels[0] is None or rels[0] != rels[1]: return None
    r = rels[0]; k = len(names)
    num = {names[i]: r[i] for i in range(k) if r[i]}; den = {names[i]: -r[k+i] for i in range(k) if r[k+i]}
    return num, den
def classify(rel):
    if rel is None: return 'unidentified'
    used = set(rel[0]) | set(rel[1])
    if 'beta4' in used: return 'beta4'
    if 'G' in used: return 'G'
    return 'pi-type'
def task(r):
    kind = r['kind']; x0, y0 = STARTS[kind]; x0 = list(x0); y0 = list(y0)
    if r.get('shift') == 'x+1': x0 = [v+1 for v in x0[:-1]] + [x0[-1]]
    elif r.get('shift') == 'y+1': y0 = [v+1 for v in y0]
    t = time.time()
    try:
        p, q = CL.walk(x0, y0, -1, r['xs'], r['ys'], [DEPTH])[DEPTH]
        g = math.gcd(p, q); p //= g; q //= g
        d100 = r['delta'][-1]
        digits = int(0.8*(1+d100)*math.log10(abs(q))) if abs(q) > 1 else 0
        out = dict(kind=kind, xs=r['xs'], ys=r['ys'], shift=r.get('shift', 'none'), delta=r['delta'], digits=digits)
        if digits < MIN_DIGITS:
            out.update(cls='too few digits', rel=None)
        else:
            rel = moebius(p, q, PI_NAMES, digits)
            if rel is None: rel = moebius(p, q, FULL[kind], digits)
            out.update(cls=classify(rel), rel=rel)
        mp.mp.dps = 30; out['L'] = mp.nstr(mp.mpf(p)/q, 20) if q else 'inf'
    except ArithmeticError:
        out = dict(kind=kind, xs=r['xs'], ys=r['ys'], shift=r.get('shift', 'none'), cls='singular', rel=None)
    out['secs'] = time.time() - t
    return out

if __name__ == "__main__":
    say(f"b4_cmf_filter.py  workers={WORKERS} scan={SCAN_JSON} depth={DEPTH}")
    # F3 Apery
    a = [F(0), F(6)]; b = [F(1), F(5)]
    for n in range(1, 150):
        c = 34*n**3+51*n**2+27*n+5
        a.append((c*a[n] - n**3*a[n-1])/(n+1)**3); b.append((c*b[n] - n**3*b[n-1])/(n+1)**3)
    rr = a[150]/b[150]
    rel = moebius(rr.numerator, rr.denominator, FULL['z3'], 300)
    ok = rel is not None and set(rel[0]) | set(rel[1]) <= {'1', 'zeta3'} and 'zeta3' in (set(rel[0]) | set(rel[1]))
    gate("F3", ok, f"relation {rel}")
    say(f"F3 Apery's recurrence at depth 150: relation {rel}")
    S = json.load(open(SCAN_JSON))
    items = [r for r in S['C2'] + S['A'] + S['A2'] if 'delta' in r and all(isinstance(x, (int, float)) and math.isfinite(x) for x in r['delta'])]
    if LIMIT_ITEMS: items = sorted(items, key=lambda r: -r['delta'][-1])[:LIMIT_ITEMS]
    say(f"{len(items)} valid trajectories from the scan")
    res = []
    with Pool(WORKERS) as pool:
        for i, o in enumerate(pool.imap_unordered(task, items, chunksize=2)):
            res.append(o)
            if (i+1) % 50 == 0: say(f"   {i+1}/{len(items)} done")
    def find(kind, xs, ys):
        return [o for o in res if o['kind'] == kind and o['xs'] == xs and o['ys'] == ys and o['shift'] == 'none']
    mp.mp.dps = 30
    f1 = find('b4', [0, 1, 1, 1, 1], [0, 0, 0, 2]); t1 = (mp.pi-2)/(8*mp.pi-12)
    ok1 = bool(f1) and f1[0]['cls'] == 'pi-type' and abs(mp.mpf(f1[0]['L']) - t1) < mp.mpf('1e-15')
    gate("F1", ok1, f"{f1[0] if f1 else 'missing'}")
    f2 = find('cat', [1, 1, 1], [0, 2]); t2 = (mp.pi-2)/(4*mp.pi-4)
    ok2 = bool(f2) and f2[0]['cls'] == 'pi-type' and abs(mp.mpf(f2[0]['L']) - t2) < mp.mpf('1e-15')
    gate("F2", ok2, f"{f2[0] if f2 else 'missing'}")
    say(f"F1 beta(4) best trajectory: class {f1[0]['cls'] if f1 else '-'}, relation {f1[0]['rel'] if f1 else '-'}")
    say(f"F2 Catalan best trajectory: class {f2[0]['cls'] if f2 else '-'}, relation {f2[0]['rel'] if f2 else '-'}")
    for kind in ('cat', 'b4'):
        R = [o for o in res if o['kind'] == kind]
        counts = {}
        for o in R: counts[o['cls']] = counts.get(o['cls'], 0) + 1
        say(f"{kind}: classes {counts}")
        for cls in (['G'] if kind == 'cat' else ['beta4']) + ['pi-type', 'unidentified']:
            C_ = sorted((o for o in R if o['cls'] == cls), key=lambda o: -o['delta'][-1])
            for o in C_[:8]:
                say(f"   {cls:12s} x {o['xs']} y {o['ys']} {o['shift']}: delta(100) {o['delta'][-1]:+.4f}  ({o['digits']} digits)  {o['rel'] if cls != 'unidentified' else ''}")
    total = sum(FAILS.values())
    json.dump(dict(stamp=STAMP, gates=FAILS, results=res, gate_fails=total), open(OUTF, "w"), indent=1, default=str)
    say(f"gates: {FAILS}")
    say(f"GATE_FAILS = {total}")
    say(f"log saved to {LOGF}")
    say(f"json saved to {OUTF}")
