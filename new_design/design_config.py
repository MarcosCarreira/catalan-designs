"""
design_config.py -- the design certified by the scripts of this folder.

Every certification script here (cat_cert4_new.py, cat_newdesign_check.py, cat_tail_opt.py,
cat_tail_edges.py, cat_cert7_new.py) reads DES and P from this file.
  DES = (n0, b0, q0, g0, h0), per unit N; a0 = b0 + q0, L0 = n0 + a0, C0 = L0 + g0, H0 = C0 + h0.
  P   = the prime modulus of cat_cert4_new.py; it must exceed H0 (gate N3 checks this).
For the paper's design REFERENCE = (96, 15, 6, 8, 5) the outputs go next to the scripts, as in the
published record runs. For any other design they go to runs/<n0>_<b0>_<q0>_<g0>_<h0>/, so the
record outputs are never overwritten.
"""
import os

DEFAULT = (36, 6, 2, 3, 2)  # design of issue #1 (fringewidth, 15 Oct 2026)
REFERENCE = (96, 15, 6, 8, 5)
# run_design_cert.py runs several designs in turn by setting the environment variable CAT_DESIGN
# ("n0,b0,q0,g0,h0") for each script it starts; without it the scripts use DEFAULT.
_ENV = os.environ.get("CAT_DESIGN")
DES = tuple(int(x) for x in _ENV.split(",")) if _ENV else DEFAULT

# prime modulus of cat_cert4_new.py: the first prime above H0 = n0 + b0 + q0 + g0 + h0
PRIMES = {(96, 15, 6, 8, 5): 131, (36, 6, 2, 3, 2): 53, (88, 14, 5, 7, 5): 127, (52, 8, 3, 4, 3): 71}
P = PRIMES[DES]

# f_low of each design from an independent implementation, for the consistency gates (T0, C6):
#   (96, 15, 6, 8, 5): cat_design_lower.py (integral 17114);
#   the others: f_low() of cat_design_search.py, exact integrals 4825/2, 57429/4, 19997/4 (2026-10-20).
FLOW_REFS = {(96, 15, 6, 8, 5): -2.275548, (36, 6, 2, 3, 2): -2.278736,
             (88, 14, 5, 7, 5): -2.270742, (52, 8, 3, 4, 3): -2.265262}
FLOW_REF = FLOW_REFS[DES]

# expected proof scores (gate T1) of cat_tail_opt.py and cat_tail_edges.py: pre-registered for the
# record design; for other designs T1 only asks for a positive score.
T1_MIN_OPT, T1_MIN_EDGES = (0.0020, 0.0025) if DES == REFERENCE else (0.0, 0.0)

def outdir(here):
    if tuple(DES) == REFERENCE:
        return here
    d = os.path.join(here, "runs", "_".join(str(x) for x in DES))
    os.makedirs(d, exist_ok=True)
    return d
