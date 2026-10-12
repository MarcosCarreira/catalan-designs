# Verification of the preprint "Catalan's constant is irrational" (OpenAI, 24 September 2026)

Marcos Costa Santos Carreira, 7 October 2026.

This folder is the verification package sent on 7 October 2026, without the covering letter and without
the copy of the preprint (linked from the top-level README). The files are unchanged.

| Path | Content |
|---|---|
| `verification_report.pdf`, `verification_report.tex` | the report (6 pages) |
| `scripts/cat_cert4.py` | rebuilds the Section 4.4 certificate (pivots mod 101, three determinants) |
| `scripts/cat_cert7.py` | rebuilds the Section 7 certificate (root brackets, value bounds, final constants) |
| `scripts/cat_padic.py` | Lemma 3.1 for 7 <= p <= 23; counting identity of Lemma 3.3; exact per-prime valuations at N = 1, 2; size of Delta_N at N = 1, 2 |
| `scripts/cat_extend.py` | Lemma 3.1 for 29 <= p <= 47; size of Delta_N at N = 1..5; valuations at the prime 2 |
| `results/*_out.json` | results written by each script on the author's computer |
| `results/*_log_*.txt` | screen output of the runs of `cat_padic.py` and `cat_extend.py` |

## Running

Python 3 with `mpmath`. No arguments; single process.

    python3 cat_cert4.py     #  ~5 s
    python3 cat_cert7.py     #  ~1 s
    python3 cat_padic.py     # ~11 min
    python3 cat_extend.py    #  ~9 min

Each script states its pass criteria in its docstring. It prints progress lines and ends with
`GATE_FAILS = <number of failed criteria>`. It writes `<name>_out.json` and a time-stamped
`<name>_log_<date>_<time>.txt` next to itself. All four runs on the author's computer ended with
`GATE_FAILS = 0`.

## Note on versions

The runs that produced `cat_cert4_out.json`, `cat_cert7_out.json` and `cat_padic_out.json` were made
before the scripts were changed to also write their screen output to a log file. That change touched only
the logging, the location of the output files and a timing note in the docstring of `cat_padic.py`. Every
computation and every pass criterion is unchanged. For `cat_padic.py` the screen output of that run is
included as `results/cat_padic_log_20261006_record.txt`. `cat_extend.py` was run in its present form.

The coefficient tables used by `cat_cert7.py`, and the pivot table used by `cat_cert4.py`, were
transcribed by hand from the preprint. The exact agreement of the rebuilt results with those tables also
confirms the transcription.
