# Catalan designs

Marcos Costa Santos Carreira, October 2026.

This repository holds the note *General designs and an irrationality exponent for Catalan's constant*
(`paper/catalan_design_extension`), its companion note *Determinant constructions for the Dirichlet beta value
at four* (`paper/beta4_determinants`), the scripts and run logs behind every computed statement in them, and the
verification report on the preprint they build on.

Archived on Zenodo: *General designs and an irrationality exponent for Catalan's constant*, doi:10.5281/zenodo.23326324
(all versions; v1 is doi:10.5281/zenodo.23326325).

The preprint is OpenAI, *Catalan's constant is irrational*, 24 September 2026:
<https://github.com/openai/math/blob/main/preprints/Catalans-constant-is-irrational-September-24-2026/paper.pdf>.
It is linked here and not copied.

## Contents of the note

1. Four conditions on the five design parameters under which the argument of the preprint proves that
   Catalan's constant G is irrational (Theorem 6, with a proof that follows the preprint section by section).
2. A second design, (96, 15, 6, 8, 5), that satisfies the four conditions with a certified gap of
   2.59e-3, about twenty times the gap of the preprint's design (Theorem 17).
3. An upper bound for the irrationality exponent of G (Section 5): the ingredients of R. K. Jha
   (doi:10.5281/zenodo.23237269), checked for the design (96, 15, 6, 8, 5), with the minors of size m < n bounded
   by an energy estimate for m nodes from the certificate of condition (C4), with sharper bounds for the row
   and column determinants (Schur polynomials; Section 5.5) and with the degree bound of S. Litsey, give
   mu(G) < 2 715 (Theorem 35); for the design of the preprint, Jha obtains 10^8 and Litsey 8 127
   (github.com/slitsey/catalan-irrationality-measure). Measurements of the determinants near G indicate that the method
   gives about 2.1e3 (Observation 37), and locate the remaining loss in the treatment of the row and column
   determinants outside the energy and in the bound for the entries (Section 5.7, Conjecture 43).
4. The best constant the real-place estimate can certify for a given design, written as an equilibrium
   problem for two measures (Section 6), with numerical values for both designs.
5. Appendix D: the preprint of Z.-W. Sun (arXiv:2609.04176) also states that Catalan's constant is irrational;
   D. Wachs (arXiv:2609.22339) shows that its final estimate does not follow from its derivation, and Appendix D
   confirms this by an exact computation for B = 20, ..., 100, so the note does not rely on that preprint.

## Contents of the companion note on beta(4)

Attempts to carry the method to beta(4) = 1 - 1/3^4 + 1/5^4 - ...: logarithmic kernels with the measure of the
preprint, a second measure with rational moments, very-well-poised forms of Rivoal and Zudilin, determinants of
half-integer arrays, a weighted arcsine measure whose determinants have entries in Q + Q beta(4) (with proofs),
Hankel determinants of very-well-poised forms, a scan of hypergeometric recurrences, and determinants with poles
at the negative squares, as in A. Fauzan's proof that zeta(5) is irrational (Section 11), whose analogue for
beta(4) is measured. None comes close to a proof. For the half-integer determinants and the weighted measure the
two bounds are measured (up to n = 40 and n = 32): the gap between them stays several times the gap of the same
construction for G, and the excess comes from the denominators of weight four. These are observations, labelled
as such in the note. Appendix A records our checks of the proofs of A. Fauzan (zeta(5)) and L. Lai (1, zeta(3),
zeta(4)), with no error found.

## Layout

| Folder | Content |
|---|---|
| `paper/` | `catalan_design_extension.tex` and `.pdf`; `beta4_determinants.tex` and `.pdf` (companion note) |
| `new_design/` | scripts for Sections 4 and 6, with the log and output of the run cited in Appendix C |
| `exponent/` | the measurements of Section 5 (`measure_scout.py` with its logs; `structure_scouts/` for Section 5.7; see `exponent/README.md`) |
| `beta4/` | the measurement scripts of the companion note with their logs and the helper modules they import; `scouts/` holds its exploration scripts (see `beta4/README.md`) |
| `sun/` | the check of Appendix D (`sun_catalan_check.py` with its log); `scouts/` holds the exploration scripts |
| `verification/` | the verification report on the preprint and its four scripts (see `verification/README.md`) |

## Scripts of the note

All scripts are Python 3, take no arguments (except the scouts), state their pass criteria in the
docstring, print progress lines, and end with `GATE_FAILS = <number of failed criteria>`. Each writes
`<name>_out.json` and `<name>_log_<date>_<time>.txt` next to itself. Every run listed below ended with
`GATE_FAILS = 0` on the author's computer. Requirements: `mpmath`, `numpy`, `scipy`; `fz5_check.py`, `b4_fauzan.py` and the scouts `fz_*_scout.py`, `cmp_content.py` also need `python-flint`.

| Script | Statement in the note | Run time of the logged run |
|---|---|---|
| `new_design/cat_newdesign_check.py` | structural conditions, ranks modulo 131, valuations (Section 4) | 5 min |
| `new_design/cat_design_search.py` | design search (Section 4.5) | 4 min |
| `new_design/cat_cert4_new.py` | condition (C3): pivots modulo 131, Appendix A | 6 s |
| `new_design/cat_tail_opt.py` | coefficients of the 46-coefficient certificate (Remark 16) | 3 s |
| `new_design/cat_tail_edges.py` | coefficients of Proposition 15 (Appendix B) | 4 s |
| `new_design/cat_cert7_new.py` | condition (C4): certificate for both coefficient files | 24 s |
| `exponent/measure_scout.py` | Observation 37 (radius on which the determinants stay small) | 12 min to 90 min per size set |
| `exponent/structure_scouts/inv_probe.py` | Observation 40, n = 96, 120, 144 (smallest singular value; no pass criteria) | 7 min |
| `new_design/cat_equilibrium.py` | equilibrium problem (Section 6) | 42 s |
| `beta4/b4_balance.py` | companion note, Observation 8, Table 1 | 5 min |
| `beta4/b4_denominators.py` | companion note, Observation 9 (needs `b4_balance.py` in the same folder) | 30 s |
| `beta4/b4_route3.py` | companion note, Observations 17, 18 and 20 (weighted measure, denominators by prime, Hankel determinants) | 2.5 min |
| `beta4/b4_cmf_scan.py` | companion note, Section 10 (hypergeometric recurrences, very-well-poised pairs) | 17 min |
| `beta4/b4_cmf_filter.py` | companion note, Observation 23 (limits of the recurrences; needs the output of `b4_cmf_scan.py`) | 4 min |
| `beta4/fz5_check.py` | companion note, Observation 26 (Fauzan's determinants for zeta(5), n = 1, 2, 3; positive control) | 5 min |
| `beta4/b4_fauzan.py` | companion note, Observations 27-29 (determinants with poles at the negative squares; common scale) | 25 s |
| `beta4/fz5_padic_check.py` | checks of Fauzan's local bounds (3.12) and Proposition 4.3 on the exact determinants, n = 1, 2, 3 (companion note, Section 11) | 6 min |
| `beta4/audit_fauzan_lai.py` | checks of the constants and lemmas of Fauzan and Lai listed in the companion note, Appendix A (runs the scripts in `beta4/fauzan_audit/` and `beta4/lai_audit/`) | 45 s |
| `sun/sun_catalan_check.py` | Appendix D (Proposition 3.1 of Sun at B = 20; Observation in Appendix D) | 20 s |

`cat_cert7_new.py` certifies the coefficient files `cat_tail_opt_out.json` and `cat_tail_edges_out.json`
committed here; their SHA-256 values are printed in the note. Running `cat_tail_opt.py` or
`cat_tail_edges.py` again overwrites these files, and with other versions of `numpy` or `scipy` the
optimiser may return different coefficients. To check the certificate of the note, run `cat_cert7_new.py`
on the committed files.

## Licence

Code: MIT (`LICENSE`). The note and the report (`paper/`, `verification/verification_report.*`):
Creative Commons Attribution 4.0 International (`LICENSE-paper.md`).

## Acknowledgement

The note, the scripts and this repository were prepared with the assistance of Claude (Anthropic).
All production runs were made by the author on his own computer.
