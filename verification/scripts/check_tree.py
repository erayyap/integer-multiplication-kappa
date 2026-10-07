"""Explicit check of one stage invocation with laminar tree side-wire groups (small h).
(1) F2: centre matrix + side rectangles == identity on the v triples;
(2) every side wire joins orthogonal spans (=> its label chain BF < x-node < y-node < AF is nested),
    and a random sample of full label chains is checked by exact rank computation mod p;
(3) data chains through all tree levels are nested; every node span is nondegenerate."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys, random, itertools
import numpy as np
sys.path.insert(0, f'{ROOT}/scripts')
from tree_groups import families, cover, anti
import check_invocation as ci  # reuses linear algebra helpers (runs its own small check on import)
h = ci.h
pS, pT = anti(h)
fS, fT = families(h, pS), families(h, pT)
T = list(itertools.combinations(range(h), 3)); idx = {t: i for i, t in enumerate(T)}; v = len(T)
def trip(K, p): return tuple(sorted(K + (p,)))
wires = []
for K, P in fS.items():
    for L, Q in fT.items():
        if len(set(K) & set(L)) == 2: continue
        n, rects = cover(K, P, L, Q)
        for (a0, a1), (b0, b1) in rects:
            wires.append((tuple(trip(K, p) for p in P[a0:a1]), tuple(trip(L, q) for q in Q[b0:b1])))
M = np.zeros((v, v), dtype=np.int64)
for Xs, Ss in wires:
    for X in Xs:
        for S in Ss:
            assert len(set(X) & set(S)) == 1
            M[idx[S], idx[X]] ^= 1
C = np.array([[len(set(S) & set(X)) % 2 for X in T] for S in T])
assert ((M + C) % 2 == np.eye(v, dtype=np.int64)).all()
print(f"tree h={h}: side wires {len(wires)}  F2 cover OK, all wires orthogonal")
random.seed(2)
for Xs, Ss in random.sample(wires, min(60, len(wires))):
    path = [ci.BF, ci.xgate(Xs, 'g'), ci.ygate(Ss, 'g'), ci.AF]
    assert all(ci.contained(U, V) for U, V in zip(path, path[1:]))
print("sampled full label chains nested OK")
def chain(fam, t, prio):
    p = min(t, key=lambda x: prio[x]); K = tuple(x for x in t if x != p); P = fam[K]
    i = P.index(p); nodes = []; node = ci.__dict__  # dummy
    from tree_groups import tree
    tr = tree(len(P))
    while True:
        lo, hi = tr[0]; nodes.append(tuple(trip(K, q) for q in P[lo:hi]))
        if tr[1] is None: break
        tr = tr[1] if i < tr[1][0][1] else tr[2]
    return nodes[::-1]  # leaf (own) first
for t in random.sample(T, 12):
    xs = chain(fS, t, pS); ys = chain(fT, t, pT)
    xpath = [ci.tens(ci.A, ci.ind[t].reshape(1, h))] + [ci.xgate(g, 'g') for g in xs] + [ci.AF]
    ypath = [ci.tens(ci.Bv, ci.ind[t].reshape(1, h)), ci.BF] + [ci.ygate(g, 'g') for g in ys[::-1]]
    for path in (xpath, ypath):
        assert all(ci.contained(U, V) for U, V in zip(path, path[1:]))
    for g in xs + ys:
        S = ci.span_t(g); assert ci.rank((S @ ci.G1 @ S.T) % ci.p) == len(g)
print("data chains through all tree levels nested; node spans nondegenerate: OK")
