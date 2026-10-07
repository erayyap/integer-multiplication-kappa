# Verification for result 001

Supporting code and data for [*Integer multiplication below n log n: improving κ from 2^-182 to 2^-59.335*](../preprints/Improving-kappa-in-integer-multiplication-below-n-log-n-October-7-2026/paper.pdf).

## Quick start

```sh
cd verification
python3 -m venv .venv && . .venv/bin/activate
pip install mpmath networkx numpy        # enough for the verifier and network checks
make verify                              # builds csrc/ helpers, reruns verify.py, diffs against data/verify_all.txt
```

`make verify` should end with `verify_all.txt reproduced exactly`. Package K's block reads:

```
== package K: bit h=(51, 50, 51) complex h=25  nu=11 log_m B<=10.2796
   1-tau = 2^-28.004  1-sigma = 2^-31.154  c = 2^-28.004  1-beta = 2^-23.26  eps = 0.199992  delta = 9.54e-07
   ...
   all constraints hold
   min g = 2^-58.329;  kappa = min g/2 (paper convention) -> kappa = 2^-59.335; any kappa < 2^-58.329 is valid
```

## Gate-level network checks

```sh
python3 scripts/check_cross.py 10 copy_sorted                    # package J (one global order), stage 1
ORDER=rot python3 scripts/check_cross.py 12 copy_sorted          # package K (cyclic orders), stage 1
ORDER=rot DUAL=1 python3 scripts/check_cross.py 12 copy_sorted   # package K, stage 2 (dual network)
python3 scripts/check_copy.py 12 copy_sorted                     # package I (per-kernel chains)
python3 scripts/check_copy.py 8 sorted_supp                      # relaxed-support labels: passes
python3 scripts/check_copy.py 8 bad_supp                         # negative control: must fail (AssertionError)
```

Each check builds a full invocation over Q^h mod 1000003 and verifies the F2 map y += x with all scratch restored, nested label chains (only centre edges decrease), nondegenerate gate labels, and the wire-count formula. Logs of the runs cited in the paper are in `data/check_*.log`.

## Layout

| Path | Contents |
|---|---|
| `scripts/verify.py` | exact constraint verifier for packages A3–K (Fractions + `mpmath.iv` interval logs) |
| `scripts/networks.py` | paper's motif-network counts (reproduces h = 100 exactly) and staged variants |
| `scripts/copy_sorted.py` | sorted-pattern rectangle cover and per-kernel chain widths (package I) |
| `scripts/cross_width.py`, `scripts/cross_var.py`, `scripts/cross_rot.py`, `scripts/rot_gaps.py` | cross-kernel chain widths, global (J) and cyclic (K) orders, width-by-gap anatomy |
| `scripts/check_copy.py`, `scripts/check_cross.py` | gate-level checkers (stage 1 and `DUAL=1`) |
| `scripts/headroom.py` | κ as a function of side wires per triple (the ≈ 2^-50.5 ceiling) |
| `scripts/sa_orders.py`, `scripts/offset_search.py`, `scripts/merge_rects.py` | searches reported as negative results |
| `scripts/groupings.py`, `scripts/tree_groups.py`, `csrc/` | grouped side gates (packages F, G) and C counters |
| `cover/`, `center/`, `center2/`, `relax/` | CP-SAT / ILP cover and centre searches, rank lower bounds |
| `design2/` | alternative per-coordinate designs: evaluator, subspace-label and complex-minrank searches |
| `data/` | cover counts (`sorted_n.txt`, `copy_n.txt`, `tree_h.txt`), checker logs, search logs, `verify_all.txt` |
| `research-notes.md` | the running research log, unedited |

Scripts locate the repository through their own path, so they can be run from any working directory. Some search scripts (`design2/`, `center/`, `cover/`) need the optional packages in `requirements.txt`.
