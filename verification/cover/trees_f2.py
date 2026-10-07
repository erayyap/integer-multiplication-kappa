"""Hierarchy (tree) product covers: tree from a complete flag of F2^L (sequence of functionals).
Tree T covers ordered cell (P,Q) iff lca(P), lca(Q) are disjoint (siblings-split at some node):
i.e. at the first functional (in flag order) non-constant on the quartet, the split is P|Q.
Find a set of trees with odd coverage for every (quartet, matching); cost ~ sum over trees of
#elements agreeing on first functional (they pay one chain per side)."""
import sys, itertools
from ortools.sat.python import cp_model
L = int(sys.argv[1]); N = 1 << L
pts = list(range(N))
def dot(l, x): return bin(l & x).count('1') & 1
# complete flags: sequences of functionals l1..lL linearly independent; tree depends on flag of spans,
# split at depth j by l_{j+1} modulo span(l1..lj) -> canonical: enumerate sequences, dedupe by node partition
def span(ls):
    s = {0}
    for l in ls: s |= {x ^ l for x in s}
    return s
seqs = []
def rec(ls):
    if len(ls) == L: seqs.append(tuple(ls)); return
    sp = span(ls)
    for l in range(1, N):
        if l not in sp: rec(ls + [l])
rec([])
def key(ls):  # tree determined by the chain of subspaces span(l1..lj)
    return tuple(frozenset(span(ls[:j])) for j in range(1, L + 1))
trees = {}
for s in seqs: trees.setdefault(key(s), s)
trees = list(trees.values())
print('flags', len(trees))
quartets = list(itertools.combinations(pts, 4))
def matchings(q):
    a, b, c, d = q
    return [((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))]
def covered(ls, q):
    for l in ls:
        v = [dot(l, x) for x in q]
        if len(set(v)) > 1:
            if sum(v) != 2: return None
            g0 = tuple(x for x, t in zip(q, v) if t == 0)
            return frozenset([g0, tuple(x for x in q if x not in g0)])
    return None
M = cp_model.CpModel(); xs = [M.NewBoolVar('') for _ in trees]
cov = {}
for ti, ls in enumerate(trees):
    for q in quartets:
        m = covered(ls, q)
        if m: cov.setdefault(m, []).append(xs[ti])
for q in quartets:
    for P, Q in matchings(q):
        m = frozenset([P, Q]); M.AddBoolXOr(cov.get(m, []))
cost = []
for ti, ls in enumerate(trees):
    pay = sum(1 for P in itertools.combinations(pts, 2) if dot(ls[0], P[0]) == dot(ls[0], P[1]))
    cost.append(pay)
M.Minimize(sum(c * x for c, x in zip(cost, xs)))
s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = float(sys.argv[2]) if len(sys.argv) > 2 else 60; s.parameters.num_workers = 8
st = s.Solve(M)
print(s.StatusName(st))
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    used = [trees[i] for i in range(len(trees)) if s.Value(xs[i])]
    print('trees', len(used), 'cost', s.ObjectiveValue(), 'C(N,2)', N * (N - 1) // 2, 'per elem', s.ObjectiveValue() * 2 / (N * (N - 1) // 2))
    for u in used: print(u)
