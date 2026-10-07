"""Scalar (m=1) labelings of small XOR circuits: is cost < W ever possible?
Circuit: list of CNOT gates (c,t) on W wires (each gate = one vertex shared by 2 wires).
Linear map must be a permutation rho.  Labels phi(v) in Q, contract equal edges.
Constraint phi(out(rho w)) - phi(in w) = 1.  Min #edges with differing values."""
import itertools, sys
import numpy as np

def lin_map(W, gates):
    M = np.eye(W, dtype=np.uint8)
    for c, t in gates:
        M[t] ^= M[c]
    return M

def perm_of(M):
    W = len(M)
    rho = {}
    for o in range(W):
        nz = np.nonzero(M[o])[0]
        if len(nz) != 1: return None
        rho[nz[0]] = o
    return rho if len(rho) == W else None

def edges(W, gates):
    # vertices: ('in',w), ('g',k), ('out',w)
    last = {w: ('in', w) for w in range(W)}
    E = []
    for k, (c, t) in enumerate(gates):
        for w in (c, t):
            E.append((last[w], ('g', k))); last[w] = ('g', k)
    for w in range(W):
        E.append((last[w], ('out', w)))
    return E

class UF:
    def __init__(s): s.p = {}; s.d = {}
    def find(s, x):
        if x not in s.p: s.p[x] = x; s.d[x] = 0; return x, 0
        if s.p[x] == x: return x, 0
        r, dd = s.find(s.p[x]); s.p[x] = r; s.d[x] += dd; return r, s.d[x]
    def union(s, a, b, w):  # phi(b) - phi(a) = w
        ra, da = s.find(a); rb, db = s.find(b)
        if ra == rb: return da + w == db
        s.p[rb] = ra; s.d[rb] = da + w - db; return True

def feasible(W, E, keep, rho):
    uf = UF()
    for i in keep:
        a, b = E[i]
        if not uf.union(a, b, 0): return False
    for w in range(W):
        if not uf.union(('in', w), ('out', rho[w]), 1): return False
    return True

def mincost(W, gates, rho):
    E = edges(W, gates)
    n = len(E)
    for k in range(n, -1, -1):  # number of equal edges
        for keep in itertools.combinations(range(n), k):
            if feasible(W, E, keep, rho):
                return n - k, E
    return None

if __name__ == "__main__":
    W = int(sys.argv[1]); G = int(sys.argv[2])
    pairs = [(c, t) for c in range(W) for t in range(W) if c != t]
    best = {}
    for gates in itertools.product(pairs, repeat=G):
        rho = perm_of(lin_map(W, gates))
        if rho is None or all(rho[w] == w for w in range(W)): continue
        r = mincost(W, gates, rho)
        if r is None: continue
        cost = r[0]
        if cost < W:
            print("SAVING", gates, rho, cost); sys.exit()
        key = tuple(sorted(rho.items()))
        best[key] = min(best.get(key, 99), cost)
    print(best)
