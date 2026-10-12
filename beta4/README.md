# Companion note: determinant constructions for beta(4)

## Measurement scripts (logged runs)

| Script | Statement | Logged run |
|---|---|---|
| `b4_balance.py` | Observation 8 and Table 1: size, denominator exponent and gap of the half-integer determinants at weights four and two, five designs, 6 <= n <= 40 | `_20261007_202313`, 5 min, `GATE_FAILS = 0` |
| `b4_denominators.py` | Observation 9: denominator exponents with and without the cancellation of G and pi^3 | `_20261007_203834`, 30 s, `GATE_FAILS = 0` |

| `b4_route3.py` | Observations 17, 18 and 20: gaps of the weighted-measure determinants (#7) and of the Catalan determinants with the same rows (#7c) up to n = 32, the preprint's determinant at its ratios, denominators by prime at n = 16, Hankel determinants of Rivoal-Zudilin moments | `_20261008_110012`, 2.5 min, `GATE_FAILS = 0` |
| `b4_cmf_scan.py` | Section 10: exponents of the pFq recurrences for Catalan's constant and beta(4), with Apery's recurrence as control; Rivoal-Zudilin pairs | `_20261008_111915`, 17 min, `GATE_FAILS = 0` |
| `b4_cmf_filter.py` | Observation 23: the limits of the trajectories of the scan, identified by PSLQ (reads `b4_cmf_scan_out_20261008_111915.json`) | `_20261008_115501`, 4 min, `GATE_FAILS = 0` |
| `fz5_check.py` | Observation 26: A. Fauzan's Hankel determinants for zeta(5) at n = 1, 2, 3 (K = 40n), primitive integer polynomial Q, log Q(zeta(5)); positive control for the measure Lambda | `_20261008_150759`, 5 min, `GATE_FAILS = 0` |
| `b4_fauzan.py` | Observations 27-29: the functionals (Z) and (B) of Definition 24 for zeta(3, 5, 7) and beta(2, 4, 6), pole steps q = 1, 2, 3, design grids, and the arcsine-family determinants on the same scale | `_20261008_153547`, 25 s, `GATE_FAILS = 0` |
| `fz5_padic_check.py` | Section 11: Fauzan's bounds (3.12) and Proposition 4.3 at every prime, on the exact determinants for n = 1, 2, 3 (needs `fz5_check.py` in the same folder) | `_20261008_160037`, 6 min, `GATE_FAILS = 0` |
| `audit_fauzan_lai.py` | Section 11: runs the checks in `fauzan_audit/` (Lemmas 3.1-3.2, 6.1, Appendices A-B, Sections 5-6 of Fauzan) and `lai_audit/` (Lemmas 2.2, 4.1, 4.4, Sections 5-6 of Lai) with 12 gates; the log contains the full output of each | `_20261008_163817`, 45 s, `GATE_FAILS = 0` |

All take no arguments and use 2 worker processes. Requirement: `mpmath`; `fz5_check.py` and `b4_fauzan.py` also need `python-flint`.

`b4_route3.py` and `b4_cmf_scan.py` import code from helper modules that must sit in the same folder:
`nu_arrays.py`, `scout_b4l.py`, `scout_b4n.py`, `full_cat.py`, `rz_scout.py`, `rz_hankel.py`, `rz_h1.py`, `b4_balance.py`
and `cmf_lite.py`; `b4_fauzan.py` imports `full_cat.py`, `scout_b4n.py`, `scout_b4l.py` and `nu_arrays.py`. These modules are kept in `beta4/` for that reason; identical copies are in `scouts/`, so that
each folder runs on its own. `cmf_lite.py` is a dependency-free implementation of the pFq conservative matrix
field of `ramanujantools` 0.0.1; gate C0 of `b4_cmf_scan.py` checks it against values computed with that package.

## Checks of Fauzan and Lai (`fauzan_audit/`, `lai_audit/`)

The scripts in these folders are run by `audit_fauzan_lai.py`; each can also be run alone from inside its folder.
`fauzan_audit/padic.py` (an earlier version of `fz5_padic_check.py`, n = 1, argument `n`) needs `fz5_check.py` copied next to it.
`lai_audit/lai_padic.py`, `lai_sum.py`, `lai_cont.py` and `lai_sec5.py` take comma-separated values of n as argument; `lai_trend.py` takes `den` or `size` followed by values of n.

## Exploration scripts (`scouts/`)

These scripts were written and run during the exploration that led to the companion note. They were run
in a sandbox, not on the author's computer, and their screen output was recorded in the author's notes;
there are no log files for them. The statements of Sections 2 to 4 and the base values of Section 6 are identifications by the PSLQ
algorithm made with these scripts, and are labelled as observations in the note. Some scripts take
arguments, shown below; the others take none.

| Script | Section of the companion note | Content |
|---|---|---|
| `scout_b4.py` | 2 | starting values of the arcsine family at weight four; PSLQ against weight-four bases |
| `scout_b4b.py` | 2 | the same with mixed-weight bases |
| `scout_b4c.py` | 2 | the sums A_n, B_n, C_n with the weights 4^k / binomial(2k, k) |
| `scout_b4d.py` | 3 | exact recurrence for the integrals J_n of theta^2 sin^n(theta) over (0, pi/2), in Q + Q pi, and one starting value |
| `scout_b4f.py` | 3 | starting values of both measures; argument `A012` or `B01` (family, then the values of r) |
| `b4g.py` | 3 | starting values of the second measure at higher precision; arguments `B<r><kind> [digits]`, e.g. `B1S1 160` |
| `rz_scout.py` | 4 | very-well-poised forms of Rivoal and Zudilin, q = 6: exact partial fractions, coefficients of beta(2) and beta(4) |
| `rz_h1.py` | 4 | elimination of beta(2) between two parameter families; denominators against size; optional argument: bound on n, exclusive (default 7) |
| `rz_h1b.py` | 4 | the same at 200 digits |
| `h2_scout.py` | 5 | first version of the half-integer determinants; arguments `ra rh rb rq n1,n2,... [W]` |
| `h2_check.py` | 5 | positive control: entries of `h2_scout.py` against the directly summed series |
| `scout_b4i.py` | 9 | starting values for odd weights in theta; arguments `w which digits`, w = 1 or 3, which = S0, S1 or P1 |
| `scout_b4j.py` | 6 | moments and starting values of the weight abs(sin theta) (theta^2 - pi theta); arguments `which digits`, which = mom, S0, S1 or P1 |
| `scout_b4k.py` | 6 | boundary sequences of the weighted measure, PSLQ; arguments `digits imax` |
| `nu_arrays.py` | 6 | exact arrays of the weighted measure and of the arcsine measure (module) |
| `nu_check.py` | 6 | positive control: `nu_arrays.py` against quadrature |
| `scout_b4l.py` | 6 | weighted-measure determinant with rows above the columns; arguments `n1,n2,... h q da b` |
| `scout_b4m.py` | 6 | arcsine array with the logarithmic kernel, PSLQ; argument `digits` |
| `scout_b4n.py` | 6 | weighted-measure determinant with overlapping rows and columns (#7) and the Catalan analogue (#7c); arguments `n1,n2,... h q a b [check]` |
| `full_cat.py` | 7 | the preprint's determinant for G with the same normalisation; arguments `n,b,q,g,h ...` or `ratio n1,n2,...` |
| `audit_r3.py` | 6.5, 6.6 | audit of the base values and of whole entries; denominators split by prime ranges; arguments `A1`, `A2`, `D` |
| `p1_profile.py` | 7 | denominator valuations binned by p/n; arguments `name n1,n2,...`, name = B7, C7, PRE, H4 or H2 (`H4` and `H2` read `b4_balance.py` from the working directory: run them from `beta4/` as `python3 scouts/p1_profile.py H4 10`) |
| `p3_small.py` | 7 | small primes: valuation of the determinant denominator against the row bound; argument `n` |
| `rz_hankel.py` | 8 | Hankel determinants of Rivoal-Zudilin moments; arguments `eta N m [check]` |
| `fz_zeta_scout.py` | 11 | Fauzan-type determinant for zeta(s); arguments `K N h r s [q] [check]` (q: pole step) |
| `fz_beta_scout.py` | 11 | its analogue for beta(s), s even; arguments `K N h r s [check]` |
| `cmp_content.py` | 11 | the arcsine-family determinants (OAI, C7, B7) with the normalization Lambda; arguments `names n1,n2,...` |

`rz_h1.py` and `rz_h1b.py` read code from `rz_scout.py`, and `h2_check.py` from `h2_scout.py`; run them
from inside `scouts/`. Several of the newer scouts import `nu_arrays.py`, `scout_b4l.py` or `scout_b4n.py` from the
same folder. Requirement: `mpmath`.
