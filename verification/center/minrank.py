"""Min F2 rank c of C (v x v, triples of [h]) with C_SS=1, C_ST=0 if |S cap T| in {0,2},
free if |S cap T|=1.  Centre loss = c*h (each centre wire loses P(x)F, dim h).
C = Y X^T: x_T, y_S in F2^c.  CP-SAT feasibility for given c."""
import sys, itertools
from ortools.sat.python import cp_model
h, c, tl = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]) if len(sys.argv) > 3 else 600
T = list(itertools.combinations(range(h), 3))
m = cp_model.CpModel()
X = [[m.NewBoolVar(f'x{s}_{r}') for r in range(c)] for s in range(len(T))]
Y = [[m.NewBoolVar(f'y{s}_{r}') for r in range(c)] for s in range(len(T))]
def prod(a, b):
    p = m.NewBoolVar(''); m.AddMultiplicationEquality(p, [a, b]); return p
# symmetry: fix first triple
nb = 0
for i, S in enumerate(T):
    for j, U in enumerate(T):
        k = len(set(S) & set(U))
        if k == 1: continue
        lits = [prod(Y[i][r], X[j][r]) for r in range(c)]
        if k == 3: m.AddBoolXOr(lits)
        else: m.AddBoolXOr(lits + [m.NewConstant(1)])
        nb += 1
sol = cp_model.CpSolver(); sol.parameters.max_time_in_seconds = tl; sol.parameters.num_workers = 8
st = sol.Solve(m)
print('h', h, 'c', c, sol.StatusName(st), flush=True)
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    ones = sum(1 for i, S in enumerate(T) for j, U in enumerate(T) if len(set(S) & set(U)) == 1
               and sum(sol.Value(Y[i][r]) * sol.Value(X[j][r]) for r in range(c)) % 2)
    tot = sum(1 for S in T for U in T if len(set(S) & set(U)) == 1)
    print('|cap|=1 ones', ones, 'of', tot)
    import json; json.dump({str(S): [[sol.Value(X[i][r]) for r in range(c)], [sol.Value(Y[i][r]) for r in range(c)]] for i, S in enumerate(T)}, open(f"data/minrank_sol_{h}_{c}.json", "w"))
    for i, S in enumerate(T[:0]):
        print(S, ''.join(str(sol.Value(X[i][r])) for r in range(c)), ''.join(str(sol.Value(Y[i][r])) for r in range(c)))
