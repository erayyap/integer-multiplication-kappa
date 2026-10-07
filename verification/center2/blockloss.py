"""Min total centre LOSS (not rank): centres = symmetric rectangles G x G, G = {S : Z <= S <= Y}.
loss(G) = rank_Q Phi[G,G], Phi_ST = |S cap T| - 1.  Need sum_k 1[S,T in G_k] = 1 (mod 2) on diag,
= 0 (mod 2) on |S cap T| in {0,2}; free on |S cap T| = 1 (side wires).  No nesting constraints."""
import itertools, sys, sympy
from ortools.sat.python import cp_model
h = int(sys.argv[1]); maxZ = int(sys.argv[2]) if len(sys.argv) > 2 else 2
pts = range(h)
T3 = list(itertools.combinations(pts, 3)); idx = {S: i for i, S in enumerate(T3)}
def rank(G):
    import numpy as np
    return int(np.linalg.matrix_rank(np.array([[len(set(a) & set(b)) - 1 for b in G] for a in G], dtype=float)))
def rank_old(G):
    M = sympy.Matrix([[len(set(a) & set(b)) - 1 for b in G] for a in G]); return M.rank()
cands = {}
for y in range(3, h + 1):
    for Y in itertools.combinations(pts, y):
        for z in range(0, min(maxZ, 3) + 1):
            for Z in itertools.combinations(Y, z):
                G = tuple(S for S in T3 if set(Z) <= set(S) <= set(Y))
                if G and G not in cands:
                    cands[G] = None
print('candidates', len(cands))
rk = {}
for G in cands:
    key = None
    rk[G] = rank(G) if len(G) < 60 else rank(G)
m = cp_model.CpModel(); x = {G: m.NewBoolVar('') for G in cands}
cover = {}
for G in cands:
    for a in G:
        for b in G:
            if idx[a] <= idx[b]: cover.setdefault((idx[a], idx[b]), []).append(x[G])
for i, S in enumerate(T3):
    for j in range(i, len(T3)):
        k = len(set(S) & set(T3[j]))
        if k == 1: continue
        lits = cover.get((i, j), [])
        par = 1 if i == j else 0
        if not lits:
            assert par == 0; continue
        m.AddBoolXOr(lits + ([m.NewConstant(1)] if par == 0 else []))
m.Minimize(sum(rk[G] * x[G] for G in cands))
s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = float(sys.argv[3]) if len(sys.argv) > 3 else 300
s.parameters.num_workers = 16
r = s.Solve(m)
print('h', h, 'status', s.StatusName(r), 'loss', s.ObjectiveValue(), 'bound', s.BestObjectiveBound(), 'paper', h * h)
for G in cands:
    if s.Value(x[G]): print(rk[G], len(G), sorted(set().union(*map(set, G))), G[:3])
