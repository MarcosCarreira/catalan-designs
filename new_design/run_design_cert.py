"""
run_design_cert.py -- runs the certification chain for each design in QUEUE, in order:
  1. cat_cert4_new.py        nonvanishing of the three determinants modulo P      (about 1 s)
  2. cat_newdesign_check.py  structural checks S, Z, B, V                         (about 2 s)
  3. cat_tail_opt.py         trial sequences, 46 coefficients                     (about 2 s)
  4. cat_tail_edges.py       trial sequences, bases at the equilibrium edges     (about 3 s)
  5. cat_cert7_new.py        rigorous certificate (reads 3 and 4) and exponent    (about 20 s)
Each script's output is shown as it runs. For each design the chain stops at the first script that
exits with an error or reports GATE_FAILS other than 0, and the next design starts. The design is
passed to the scripts through the environment variable CAT_DESIGN (read by design_config.py).
At the end: a summary line per script and design, the exponent lines (E1) of cat_cert7_new.py, and
a table of the certified gap and exponent bound per design (edge-tail input).

RUN: python3 run_design_cert.py   (no arguments; needs numpy, scipy, mpmath; single process each)
TIME: about 30-60 s per design, plus the V check of cat_newdesign_check.py (about 5 min at n0 = 88).
      For the present QUEUE: about 7 min.
OUTPUT: screen + run_design_cert_log_<YYYYMMDD_HHMMSS>.txt in runs/; each script writes its own
        log and JSON in design_config.outdir() (runs/<design>/ for any design other than
        (96, 15, 6, 8, 5)).
"""
import os, sys, subprocess, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
QUEUE = [(36, 6, 2, 3, 2), (88, 14, 5, 7, 5), (52, 8, 3, 4, 3)]   # top three of cat_mu_search.py (21 Oct 2026)
OUT = os.path.join(HERE, "runs")
os.makedirs(OUT, exist_ok=True)

CHAIN = ["cat_cert4_new.py", "cat_newdesign_check.py", "cat_tail_opt.py",
         "cat_tail_edges.py", "cat_cert7_new.py"]

T0 = time.time()
STAMP = time.strftime('%Y%m%d_%H%M%S')
LOGPATH = os.path.join(OUT, 'run_design_cert_log_' + STAMP + '.txt')
LOGF = open(LOGPATH, 'w')
def log(s):
    line = f"[{time.time()-T0:6.1f}s] {s}"
    print(line, flush=True); LOGF.write(line + '\n'); LOGF.flush()

import re
summary, expo, table = [], [], []
for des in QUEUE:
    tag = ",".join(str(x) for x in des)
    env = dict(os.environ, CAT_DESIGN=tag)
    log(f"##### design {des}")
    ok_all = True
    for i, name in enumerate(CHAIN, 1):
        log(f"=== {des} {i}/{len(CHAIN)} {name}")
        t = time.time()
        proc = subprocess.Popen([sys.executable, os.path.join(HERE, name)], cwd=HERE, env=env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        fails = None
        for line in proc.stdout:
            print("    " + line, end="", flush=True)
            if line.startswith("GATE_FAILS = "): fails = int(line.split("=")[1])
            if name == "cat_cert7_new.py" and ("E1 exponent" in line or "rho_1 <=" in line or "E2 " in line):
                expo.append(f"{des}: " + line.rstrip().split("] ", 1)[-1])
                m = re.search(r"gamma >= ([0-9.e-]+), mu <= ([0-9.]+)", line)
                if m and expo[-2:] and "cat_tail_edges_out.json" in expo[-2]:
                    table.append((des, m.group(1), m.group(2)))
        rc = proc.wait()
        ok = rc == 0 and fails == 0
        summary.append(f"{str(des):20s} {name:24s} exit {rc}, GATE_FAILS {fails}, {time.time()-t:5.1f}s {'ok' if ok else 'FAILED'}")
        LOGF.write(f"{des} {name}: exit {rc}, GATE_FAILS {fails}\n")
        if not ok:
            ok_all = False
            log(f"STOP: {name} did not pass for {des}; later scripts for this design not run")
            break
    log(f"design {des}: chain {'PASSED' if ok_all else 'FAILED'}")

log("SUMMARY")
for s_ in summary: log("  " + s_)
for s_ in expo: log("  " + s_)
log("CERTIFIED (edge-tail input): design, gap gamma >=, mu <=")
for des, g, mu in table: log(f"  {str(des):20s} {g:>16s} {mu:>12s}")
LOGF.close()
print(f"log saved to {LOGPATH}")
