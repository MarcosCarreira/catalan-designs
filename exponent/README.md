# Section 5: an irrationality exponent

| Script | Statement | Logged runs |
|---|---|---|
| `measure_scout.py` | Observation 37: the exact polynomial Delta_N(X), its expansion at G, and the largest radius d* on which abs(Delta_N(x)) <= 2 abs(Delta_N(G)); prints rho_n = -log(d*)/n | `_20261008_183830` (ratio 48:7:4:4:2, n = 24 to 120, 735 s); `_20261008_192255` (ratio 96:15:6:8:5, n = 48 to 144, 1969 s); `_20261008_220701` (ratio 96:15:6:8:5, n = 192, 5433 s) |

Settings (`RATIO`, `SIZES`) are constants at the top of `measure_scout.py`; the logged runs used `RATIO = OURS` with
`SIZES = [48, 72, 96, 120, 144]` and `[192]`, and the earlier version of the script with the ratio 48:7:4:4:2 and
sizes 24 to 120. No arguments. `full_cat.py` (the construction of the matrices, also in `beta4/`) must be in the same
folder. Requirements: `python-flint`, `mpmath`.

## Structure near the constant (Section 5.7)

`structure_scouts/` holds the scripts behind Observations 40 and 42. They measure without pass criteria.

| Script | Statement | Runs |
|---|---|---|
| `inv_probe.py` | Observation 40: log of the Frobenius norm of A^-1 in ball arithmetic, an upper bound for -log(sigma_n)/n within log(n)/n; no arguments, designs set in `DESIGNS` | `inv_probe_log_20261009_105215.txt` (n = 96, 120, 144, 397 s); n = 24, 48, 72 in the sandbox |
| `sv_probe.py` | Observation 40: all singular values of A and their gaps; argument `n,b,q,g,h` and optional digits | sandbox, n = 24, 48 |
| `gcol_probe.py` | Observation 40: u^T B v on the bottom singular vectors; arguments `n,b,q,g,h,digits` | sandbox, n = 24, 48 |
| `lam_probe.py` | Observation 40: spacing of the zeros of Delta_N(X) from the expansion at G (imports `measure_scout.py`) | sandbox, n = 24, 48 |
| `theta_probe2.py`, `alt_probe.py` | Observation 42: -F'(1) by the envelope formula, r_row and r_col (arguments: designs `n,b,q,g,h`; load `tp_head.py`, which reads `cat_equilibrium.py` from the working directory) | sandbox |
| `theta_probe.py` | F(theta) on a grid of theta, for the finite-difference check of Observation 42 | sandbox |

The scouts other than `inv_probe.py` import `full_cat.py`, `measure_scout.py` or `cat_equilibrium.py`
(in `exponent/` and `new_design/`) by module name: run them from a folder that holds those files, or put
those folders on `PYTHONPATH`. Requirements: `python-flint` (`inv_probe.py`), `mpmath`, `numpy`.
