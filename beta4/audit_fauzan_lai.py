"""audit_fauzan_lai.py -- runs the checks of A. Fauzan, "zeta(5) is irrational" (Zenodo, 17 Sep 2026) and of L. Lai,
arXiv:2610.06838, kept in fauzan_audit/ and lai_audit/, and applies pass criteria to their output.

Run:  python3 audit_fauzan_lai.py        (no arguments)
Needs: python-flint, mpmath, sympy, numpy; the folders fauzan_audit/ and lai_audit/ next to this file.
Writes audit_fauzan_lai_log_<stamp>.txt (with the full output of every check) next to this file.
Expected runtime: about 1.5 minutes (the checks run one after the other, one process each).
The local estimates of Fauzan on the exact determinants are checked separately by fz5_padic_check.py.

Gates
 F1 Lemma 6.1 of Fauzan: potential inequality and energy constant (lemma61.py)
 F2 Appendix B: int_3^20 R/x^3 equals the stated rational exactly (consts2.py)
 F3 Appendix B: I_out = 127751/96000 and int d = 9/640 exactly (iout.py); A* = I_out + int_3^20 - 2689/48000 exactly
 F4 Lemmas 3.1-3.2: kappa_d valuations, distribution formula on polynomials and at poles (v_p >= 43), C_p in p^5 Z_p (lemma32.py)
 F5 Section 5 tail: -1/2 <= Q <= 13/8, mean of P = 2923/240, |C| < 16, P(20) = 37/2, C(20) = 0 (tail_check.py)
 F6 Section 6: (6.8) <= 0, (6.11), (6.12), the constant 652, (6.15) (sec6.py)
 F7 (5.7): |gamma_p^in - p Gamma(K/p)| <= 5M at K = 320000, M = 40 (gamma_in.py)
 F8 outer limit: |gamma_p^out + K(R0 - d)(p/K)| <= 10 at K = 320000 (gamma_out.py)
 L1 Lai Lemma 2.2 (series = formula) and Proposition 3.1 at q = 19 (lai_check.py, first part)
 L2 Lai Lemmas 4.1 and 4.4 on the exact J_n, n = 8, 12, 16, 20, 24 (lai_padic.py)
 L3 Lai Section 5: Table 1 integral, finite checks of Lemma 5.4, Section-5 bound <= Lemma 4.4 bound at every prime (lai_sec5.py)
 L4 Lai Section 6: lambda0, c1 < 0, c2 > 0, c3 = closed form, G - lambda0 y <= c1, Lemma 6.2, identity (6.9), Lemma 6.1 (lai_sec6.py)
"""
import os, sys, re, time, datetime, subprocess
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
STAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
LOGF = os.path.join(HERE, f"audit_fauzan_lai_log_{STAMP}.txt")
T0 = time.time(); _log = open(LOGF, "w")
def say(s):
    line = f"[{time.time()-T0:7.1f}s] {s}"; print(line, flush=True); _log.write(line + "\n"); _log.flush()
FAILS = {}
def gate(name, ok, info=""):
    FAILS[name] = FAILS.get(name, 0) + (0 if ok else 1)
    say(f"{'pass' if ok else 'FAIL'} {name} {info}")
def run(folder, script, *args):
    say(f"running {folder}/{script} {' '.join(args)}")
    out = subprocess.run([sys.executable, script, *args], cwd=os.path.join(HERE, folder), capture_output=True, text=True)
    for l in (out.stdout + out.stderr).splitlines(): _log.write("    | " + l + "\n")
    _log.flush()
    if out.returncode != 0: say(f"   exit code {out.returncode}")
    return out.stdout
def nums(s): return [float(x) for x in re.findall(r'[-+]?\d+\.?\d*(?:e[-+]?\d+)?', s)]

o = run('fauzan_audit', 'lemma61.py'); gate("F1", o.count("ok: True") == 2)
o = run('fauzan_audit', 'consts2.py', '3', '20'); gate("F2", "equal exactly: True" in o)
o = run('fauzan_audit', 'iout.py')
A = Fr(127751, 96000) + Fr(322437603634266857629, 7535670527041937280000) - Fr(2689, 48000)
gate("F3", "127751/96000: True" in o and "9/640: True" in o and A == Fr(9928298118277006344769, 7535670527041937280000))
o = run('fauzan_audit', 'lemma32.py')
vmins = [int(x) for x in re.findall(r'min v_p\(difference\) = (\d+)', o)]; cps = [int(x) for x in re.findall(r'v_p\(C_p\) = (\d+)', o)]
gate("F4", "violations [] 0" in o and "(exact): True" in o and len(vmins) == 3 and min(vmins) >= 43 and len(cps) == 12 and min(cps) >= 5,
     f"pole identity v_p {vmins}, C_p valuations {cps}")
o = run('fauzan_audit', 'tail_check.py')
q = re.search(r'min ([-\d.e]+)\s+max ([-\d.e]+)', o); c = re.search(r'max \|C\| sampled ([\d.]+)', o); pm = re.search(r'mean of P .* = ([\d.]+)', o)
gate("F5", q and float(q.group(1)) >= -0.5 and float(q.group(2)) <= 13/8 and c and float(c.group(1)) < 16 and pm and abs(float(pm.group(1)) - 2923/240) < 1e-5
     and "P(20) = 37/2  C(20) = 0" in o, f"Q in [{q.group(1)}, {q.group(2)}], max|C| {c.group(1)}, mean P {pm.group(1)}")
o = run('fauzan_audit', 'sec6.py')
m68 = re.search(r'\(6\.8\): max of .*?: ([-\d.e]+)', o); m611 = re.search(r'\(6\.11\): max ratio .*?: ([\d.e-]+)', o)
gate("F6", m68 and float(m68.group(1)) <= 0 and m611 and float(m611.group(1)) <= 1 and "(6.12): violations 0" in o and "= 652.0" in o and o.count("ok=True") == 4)
o = run('fauzan_audit', 'gamma_in.py'); d = [float(x) for x in re.findall(r'diff=([-\d.]+)', o)]
gate("F7", len(d) >= 6 and max(abs(x) for x in d) <= 200, f"differences {d}")
o = run('fauzan_audit', 'gamma_out.py'); d = [float(x) for x in re.findall(r'diff=([-\d.]+)', o)]
gate("F8", len(d) >= 5 and max(abs(x) for x in d) <= 10, f"differences {d}")
o = run('lai_audit', 'lai_check.py')
r = re.search(r'series vs formula: ([\d.e-]+)', o)
gate("L1", r and float(r.group(1)) < 1e-40 and o.count("J_4 != 0: True") == 4, f"Lemma 2.2 difference {r.group(1) if r else '?'}")
o = run('lai_audit', 'lai_padic.py', '8,12,16,20,24'); gate("L2", o.count("violations [];") == 5)
o = run('lai_audit', 'lai_sec5.py', '30,60,100'); gate("L3", "-657/70: True" in o and o.count("bound: [] (0); finite-check failures in Lemma 5.4: [] (0)") == 3)
o = run('lai_audit', 'lai_sec6.py')
c3 = re.search(r'c3 = ([\d.]+)', o); c3p = re.search(r'closed form = ([\d.]+)', o)
mx = [float(x) for x in re.findall(r'\(must be <= 0\)', o) and re.findall(r': ([-\d.e]+)\s*\(must be <= 0\)', o)]
h69 = re.findall(r'Hankel det of moments = ([\d.e-]+)\s+prod k!\^2/4\^k = ([\d.e-]+)', o)
ctr = re.findall(r'contour \(([-\d.e]+)[^)]*\)\s+series ([-\d.e]+)', o)
gate("L4", c3 and c3p and abs(float(c3.group(1)) - float(c3p.group(1))) < 1e-7 and abs(float(c3.group(1)) - 31.614437) < 1e-5 and len(mx) == 2 and max(mx) <= 0
     and len(h69) == 4 and all(abs(float(a) - float(b)) < 1e-12 for a, b in h69) and len(ctr) == 3 and all(abs(float(a) - float(b)) <= 1e-12*abs(float(b)) for a, b in ctr),
     f"c3 {c3.group(1) if c3 else '?'}")
total = sum(FAILS.values())
say(f"gates: {FAILS}")
say(f"GATE_FAILS = {total}")
say(f"log saved to {LOGF}")
