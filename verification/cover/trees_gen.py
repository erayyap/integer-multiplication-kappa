"""All rooted binary trees on n labelled leaves; choose odd-cover set minimizing
cost = sum_T [ 2*#{P: lca_T(P) != root} + #useful rects ] (rect (c1,c2) & (c2,c1) at a node if both children have >=2 leaves)."""
import sys, itertools
from ortools.sat.python import cp_model
n = int(sys.argv[1]); tl = float(sys.argv[2]) if len(sys.argv) > 2 else 120
def trees(leaves):
    if len(leaves) == 1: yield leaves[0]; return
    first, rest = leaves[0], leaves[1:]
    # split: subset S containing first vs complement nonempty
    for r in range(0, len(rest)):
        for S in itertools.combinations(rest, r):
            A = (first,) + S; B = tuple(x for x in rest if x not in S)
            for ta in trees(A):
                for tb in trees(B): yield (ta, tb)
def leaves(t): return (t,) if isinstance(t, int) else leaves(t[0]) + leaves(t[1])
def splits(t):  # list of (A,B) child leaf sets at internal nodes
    if isinstance(t, int): return []
    return [(frozenset(leaves(t[0])), frozenset(leaves(t[1])))] + splits(t[0]) + splits(t[1])
T = list(trees(tuple(range(n))))
print('trees', len(T), flush=True)
pairs = list(itertools.combinations(range(n), 2))
M = cp_model.CpModel(); xs = [M.NewBoolVar('') for _ in T]
cov = {}; cost = []
for i, t in enumerate(T):
    sp = splits(t); root = sp[0]
    pay = sum(1 for P in pairs if set(P) <= root[0] or set(P) <= root[1])
    rects = sum(2 for A, B in sp if len(A) >= 2 and len(B) >= 2)
    cost.append(2 * pay + rects)
    for A, B in sp:
        for P in itertools.combinations(sorted(A), 2):
            for Q in itertools.combinations(sorted(B), 2):
                cov.setdefault(frozenset([P, Q]), []).append(xs[i])
for q in itertools.combinations(range(n), 4):
    a, b, c, d = q
    for P, Q in (((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))):
        M.AddBoolXOr(cov.get(frozenset([P, Q]), []))
M.Minimize(sum(c * x for c, x in zip(cost, xs)))
s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = 8
st = s.Solve(M); print(s.StatusName(st), 'cost', s.ObjectiveValue(), 'bound', s.BestObjectiveBound(), 'C(n,2)', len(pairs))
for i in range(len(T)):
    if s.Value(xs[i]): print(T[i], cost[i])
