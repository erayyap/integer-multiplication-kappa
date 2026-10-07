"""CP-SAT: min (#rects + sum widths) exact rect cover of disjoint-pair cells, nodes = K(A) or star(x,B)."""
import sys, itertools
from ortools.sat.python import cp_model
n = int(sys.argv[1]); tl = float(sys.argv[2]) if len(sys.argv) > 2 else 120
pairs = list(itertools.combinations(range(n), 2)); pid = {p: i for i, p in enumerate(pairs)}
nodes = {}
for r in range(2, n + 1):
    for A in itertools.combinations(range(n), r):
        mem = frozenset(pid[p] for p in itertools.combinations(A, 2)); nodes.setdefault(mem, set(A))
for x in range(n):
    rest = [y for y in range(n) if y != x]
    for r in range(1, n):
        for B in itertools.combinations(rest, r):
            mem = frozenset(pid[tuple(sorted((x, b)))] for b in B); nodes.setdefault(mem, {x, *B})
NL = list(nodes); sup = [nodes[m] for m in NL]
print('nodes', len(NL), flush=True)
rects = [(i, j) for i in range(len(NL)) for j in range(len(NL)) if not (sup[i] & sup[j])]
print('rects', len(rects), flush=True)
M = cp_model.CpModel()
x = [M.NewBoolVar('') for _ in rects]
cov = {}
for k, (i, j) in enumerate(rects):
    for a in NL[i]:
        for b in NL[j]: cov.setdefault((a, b), []).append(x[k])
for P in range(len(pairs)):
    for Q in range(len(pairs)):
        if set(pairs[P]) & set(pairs[Q]): continue
        M.AddExactlyOne(cov[(P, Q)])
us = [M.NewBoolVar('') for _ in NL]; ut = [M.NewBoolVar('') for _ in NL]
byi = {}; byj = {}
for k, (i, j) in enumerate(rects): byi.setdefault(i, []).append(x[k]); byj.setdefault(j, []).append(x[k])
for i in range(len(NL)):
    for u, by in ((us, byi), (ut, byj)):
        lits = by.get(i, [])
        for l in lits: M.AddImplication(l, u[i])
        M.AddBoolOr(lits + [u[i].Not()])
cost = [sum(x)]
W1 = True
contain = {P: [i for i in range(len(NL)) if P in NL[i]] for P in range(len(pairs))}
for P in range(len(pairs)):
    C = contain[P]
    for u in (us, ut):
        for a in range(len(C)):
            for b in range(a+1, len(C)):
                i1, j1 = C[a], C[b]
                if not (NL[i1] <= NL[j1] or NL[j1] <= NL[i1]): M.AddBoolOr([u[i1].Not(), u[j1].Not()])
cost.append(2 * len(pairs))
M.Minimize(sum(cost))
s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = 8
st = s.Solve(M)
print('n', n, s.StatusName(st), 'obj', s.ObjectiveValue(), 'bound', s.BestObjectiveBound())
used = [rects[k] for k in range(len(rects)) if s.Value(x[k])]
def show(i):
    m = NL[i]; return sorted(pairs[p] for p in m)
for i, j in used: print(show(i), 'x', show(j))
