"""Min-loss F2 centre decomposition search (small h) with CP-SAT.
Rect types (both sides contain kernel point i):
  sym : {S: i in S, S-i in Z} x same,            loss = rank Gram block
  asym: {S: i in S, S-i in A} x {T: i in T, T-i in B}
Requirement (ordered pairs): S=T odd; |S cap T| = 2 even; |S cap T|=1 free; disjoint untouched."""
import sys, itertools
from fractions import Fraction as Fr
import numpy as np
from ortools.sat.python import cp_model
h = int(sys.argv[1]); mode = sys.argv[2] if len(sys.argv) > 2 else 'sym'
tl = float(sys.argv[3]) if len(sys.argv) > 3 else 60
T = list(itertools.combinations(range(h), 3))
G1 = np.eye(h) - np.ones((h, h)) / 9
def ind(S): v = np.zeros(h); v[list(S)] = 1; return v
def grank(Ss, Ts):
    if not Ss or not Ts: return 0
    M = np.array([ind(S) for S in Ss]) @ G1 @ np.array([ind(S) for S in Ts]).T
    return np.linalg.matrix_rank(M, tol=1e-9)
rects = []
for i in range(h):
    rest = [j for j in range(h) if j != i]
    def fam(Z): return [tuple(sorted((i,) + p)) for p in itertools.combinations(sorted(Z), 2)]
    subsets = [Z for r in range(2, h) for Z in itertools.combinations(rest, r)]
    if mode == 'sym':
        for Z in subsets:
            F = fam(Z); rects.append((F, F, grank(F, F)))
    else:
        for A in subsets:
            for B in subsets:
                FA, FB = fam(A), fam(B); r = grank(FA, FB)
                if r == 0: continue
                # skip rects that touch no constrained pair
                rects.append((FA, FB, r))
print('rects', len(rects), flush=True)
cons = {}
for S in T:
    for U in T:
        k = len(set(S) & set(U))
        if k in (2, 3): cons[(S, U)] = 1 if k == 3 else 0
m = cp_model.CpModel(); x = [m.NewBoolVar(f'x{k}') for k in range(len(rects))]
touch = {key: [] for key in cons}
for k, (FA, FB, r) in enumerate(rects):
    for S in FA:
        for U in FB:
            if (S, U) in touch: touch[(S, U)].append(x[k])
one = m.NewConstant(1)
for key, b in cons.items():
    lits = touch[key]
    if b == 1: m.AddBoolXOr(lits)
    else: m.AddBoolXOr(lits + [one])
m.Minimize(sum(r * x[k] for k, (_, _, r) in enumerate(rects)))
sol = cp_model.CpSolver(); sol.parameters.max_time_in_seconds = tl; sol.parameters.num_workers = 16
st = sol.Solve(m)
print('h', h, 'status', sol.StatusName(st), 'loss', sol.ObjectiveValue(), 'bound', sol.BestObjectiveBound(), 'paper h^2 =', h * h)
used = [rects[k] for k in range(len(rects)) if sol.Value(x[k])]
print('rects used', len(used), 'sizes', sorted((len(a), len(b), r) for a, b, r in used)[:40])
