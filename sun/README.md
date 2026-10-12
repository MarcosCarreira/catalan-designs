# Appendix D: the preprint of Z.-W. Sun

| Script | Statement | Logged run |
|---|---|---|
| `sun_catalan_check.py` | Appendix D: the identity of Proposition 3.1 of Sun checked exactly at B = 20 (gate S1), and the values of (log D_B + log abs(qhat_B(G)))/B^2 for B = 20, 40, 60, 80, 100 over all admissible row sets A (Observation in Appendix D) | `_20261008_165354`, 20 s, `GATE_FAILS = 0` |

No arguments. Requirements: `python-flint`, `mpmath`.

`scouts/` holds the sandbox scripts that preceded it: `sun_prop31.py` (Proposition 3.1 at B = 20), `sun_parts.py`
(sizes of the factors of qhat_B; argument: values of B), `sun_check.py` (the measured quantity; argument: values of B).
