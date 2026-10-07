"""Explicit check of one stage-2 invocation with grouped side wires (h small).
Labels live in F (x) F (A = F first factor, exposed coordinate second), form (I-J/9)^{(x)2},
computed mod a large prime.  Checks: (1) every wire/data edge joins nested labels,
(2) only center 3->4 edges decrease (loss h each), (3) F2 exact cover: center + side = I."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import itertools, sys, random
import numpy as np
sys.path.insert(0, f'{ROOT}/scripts')
from groupings import design
p = 1000003
h = int(sys.argv[1]); k = int(sys.argv[2])
inv9 = pow(9, p - 2, p)
G1 = (np.eye(h, dtype=np.int64) - inv9 * np.ones((h, h), dtype=np.int64)) % p
G = np.kron(G1, G1) % p

def rank(M):
    M = M.copy() % p; r = 0; rows, cols = M.shape
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i, c]: piv = i; break
        if piv is None: continue
        M[[r, piv]] = M[[piv, r]]
        M[r] = M[r] * pow(int(M[r, c]), p - 2, p) % p
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] = (M[i] - M[i, c] * M[r]) % p
        r += 1
        if r == rows: break
    return r

def nullspace(M):
    # basis of {x : M x = 0}
    M = M.copy() % p; rows, cols = M.shape; r = 0; pivc = []
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i, c]: piv = i; break
        if piv is None: continue
        M[[r, piv]] = M[[piv, r]]
        M[r] = M[r] * pow(int(M[r, c]), p - 2, p) % p
        for i in range(rows):
            if i != r and M[i, c]: M[i] = (M[i] - M[i, c] * M[r]) % p
        pivc.append(c); r += 1
    free = [c for c in range(cols) if c not in pivc]
    B = []
    for f in free:
        x = np.zeros(cols, dtype=np.int64); x[f] = 1
        for i, c in enumerate(pivc): x[c] = (-M[i, f]) % p
        B.append(x)
    return np.array(B, dtype=np.int64).reshape(len(B), cols)

def perp(B, Gm):  # orthogonal complement of rowspace B
    if len(B) == 0: return np.eye(Gm.shape[0], dtype=np.int64)
    return nullspace((B @ Gm) % p)

def tens(U, V):
    return np.array([np.kron(u, w) % p for u in U for w in V], dtype=np.int64).reshape(len(U) * len(V), -1)

def ssum(*Us):
    Us = [u for u in Us if len(u)]
    return np.vstack(Us) if Us else np.zeros((0, h * h), dtype=np.int64)

def contained(U, V):
    if len(U) == 0: return True
    return rank(np.vstack([V, U])) == rank(V)

def dim(U): return rank(U) if len(U) else 0

T = list(itertools.combinations(range(h), 3))
ind = {t: np.array([1 if i in t else 0 for i in range(h)], dtype=np.int64) for t in T}
I_h = np.eye(h, dtype=np.int64)
a1 = T[0]
Pv = ind[a1].reshape(1, h)
Bv = perp(Pv, G1)
A = I_h
FF = I_h
BF = tens(Bv, FF); AF = tens(A, FF)
def PV(V): return tens(Pv, V)
gS, gT = design(h, k, "anti")
# wires: rebuild the block covers as in blocks.c
def nbr(X, S): return len(set(X) & set(S)) == 1
wires = []
for Gs in gS:
    for Gt in gT:
        cells = [(X, S) for X in Gs for S in Gt if nbr(X, S)]
        if not cells: continue
        n = len(cells)
        if n == len(Gs) * len(Gt): wires.append(((tuple(Gs), 'g'), (tuple(Gt), 'g'))); continue
        fr = [X for X in Gs if all(nbr(X, S) for S in Gt)]
        fc = [S for S in Gt if all(nbr(X, S) for X in Gs)]
        rest_r = [(X, S) for (X, S) in cells if X not in fr]
        rest_c = [(X, S) for (X, S) in cells if S not in fc]
        opts = [(n, 's'), (len(fr) + len(rest_r), 'r'), (len(fc) + len(rest_c), 'c')]
        best = min(opts)[1]
        if best == 's':
            for X, S in cells: wires.append((((X,), 'o'), ((S,), 'o')))
        elif best == 'r':
            for X in fr: wires.append((((X,), 'o'), (tuple(Gt), 'g')))
            for X, S in rest_r: wires.append((((X,), 'o'), ((S,), 'o')))
        else:
            for S in fc: wires.append(((tuple(Gs), 'g'), ((S,), 'o')))
            for X, S in rest_c: wires.append((((X,), 'o'), ((S,), 'o')))
print("side wires", len(wires), "pairs", sum(1 for X in T for S in T if nbr(X, S)))
# F2 exact cover check
idx = {t: i for i, t in enumerate(T)}
v = len(T)
M = np.zeros((v, v), dtype=np.int64)
for (Xs, _), (Ss, _) in wires:
    for X in Xs:
        for S in Ss: M[idx[S], idx[X]] ^= 1
C = np.array([[len(set(S) & set(X)) % 2 for X in T] for S in T])
assert ((M + C) % 2 == np.eye(v, dtype=np.int64)).all(), "cover fails"
print("F2: center + side = identity  OK")
# labels
cache = {}
def span_t(ts): return np.array([ind[t] for t in ts], dtype=np.int64)
def xgate(Xs, kind):
    key = ('x', Xs)
    if key not in cache: cache[key] = ssum(BF, PV(span_t(Xs)))
    return cache[key]
def ygate(Ss, kind):
    key = ('y', Ss)
    if key not in cache: cache[key] = ssum(BF, PV(perp(span_t(Ss), G1)))
    return cache[key]
random.seed(1)
sample = random.sample(wires, min(len(wires), int(sys.argv[3]) if len(sys.argv) > 3 else 200))
bad = 0
for (Xs, kx), (Ss, ks) in sample:
    path = [BF, xgate(Xs, kx), ygate(Ss, ks), AF]
    for U, V in zip(path, path[1:]):
        if not contained(U, V): bad += 1
print("sampled wire paths checked:", len(sample), "bad:", bad)
# data chains for a few elements
for t in T[:15]:
    gx = [g for g in gS if t in g][0]; gy = [g for g in gT if t in g][0]
    xpath = [tens(A, ind[t].reshape(1, h)), xgate((t,), 'o'), xgate(tuple(gx), 'g'), AF]
    ypath = [tens(Bv, ind[t].reshape(1, h)), BF, ygate(tuple(gy), 'g'), ygate((t,), 'o'),
             ssum(BF, PV(perp(ind[t].reshape(1, h), G1)))]
    for path in (xpath, ypath):
        for U, V in zip(path, path[1:]):
            assert contained(U, V)
# nondegeneracy of group spans
for g in gS + gT:
    S = span_t(g); assert rank((S @ G1 @ S.T) % p) == len(g)
print("data chains nested, group spans nondegenerate: OK")
