"""Cross-kernel copy chains for the sorted cover.

A copy wire (source side) only needs its gate labels span(t_G) to be nested in time order;
nothing forces all nodes of one chain to be in the same kernel family.  For each triple S
take all its source nodes over its 3 kernels, order them by span inclusion over Q, and
compute the Dilworth width of that poset (min number of chains).  Same for targets
(labels span(t_H)^perp, inclusion reversed -> same poset on spans).
Compare with the per-kernel widths used by package I.
"""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys, itertools
from collections import defaultdict
import networkx as nx
sys.path.insert(0, f'{ROOT}/scripts')
import copy_sorted as cs

p = 1000003


def rref(rows, h):
    """basis (list of (pivot, row)) of span of rows mod p, fully reduced."""
    basis = []
    for r in rows:
        r = [x % p for x in r]
        for piv, b in basis:
            if r[piv]:
                f = r[piv]; r = [(x - f * y) % p for x, y in zip(r, b)]
        nz = next((i for i, x in enumerate(r) if x), None)
        if nz is None: continue
        inv = pow(r[nz], p - 2, p); r = [(x * inv) % p for x in r]
        nb = []
        for piv, b in basis:
            if b[nz]:
                f = b[nz]; b = [(x - f * y) % p for x, y in zip(b, r)]
            nb.append((piv, b))
        basis = nb + [(nz, r)]
    return basis


def inspan(vec, basis):
    r = list(vec)
    for piv, b in basis:
        if r[piv]:
            f = r[piv]; r = [(x - f * y) % p for x, y in zip(r, b)]
    return not any(r)


def dilworth(items, leq):
    L = len(items)
    Gb = nx.Graph(); Gb.add_nodes_from([('l', i) for i in range(L)] + [('r', i) for i in range(L)])
    for i in range(L):
        for j in range(L):
            if i != j and leq(items[i], items[j]) and (not leq(items[j], items[i]) or i < j):
                Gb.add_edge(('l', i), ('r', j))
    M = nx.bipartite.hopcroft_karp_matching(Gb, top_nodes=[('l', i) for i in range(L)])
    return L - len(M) // 2


def run(h):
    n = h - 1
    st = cs.structure(n)
    mem = st['members']
    srcn = {g for g, _ in st['rects']}; tgtn = {hh for _, hh in st['rects']}
    ind = lambda t: [1 if i in t else 0 for i in range(h)]
    trip_nodes = {}  # (q, nd) -> frozenset of triples
    for q in range(h):
        others = [i for i in range(h) if i != q]
        for nd in srcn | tgtn:
            trip_nodes[(q, nd)] = frozenset(tuple(sorted((q, others[a], others[b]))) for a, b in mem[nd])
    basis = {k: rref([ind(t) for t in v], h) for k, v in trip_nodes.items()}
    cache = {}
    def leq(a, b):
        key = (a, b)
        if key not in cache:
            cache[key] = trip_nodes[a] <= trip_nodes[b] or all(inspan(ind(t), basis[b]) for t in trip_nodes[a])
        return cache[key]
    by_src = defaultdict(set); by_tgt = defaultdict(set)
    for q in range(h):
        for nd in srcn:
            for t in trip_nodes[(q, nd)]: by_src[t].add((q, nd))
        for nd in tgtn:
            for t in trip_nodes[(q, nd)]: by_tgt[t].add((q, nd))
    per_kernel = 2 * h * (st_total := sum(len(c) for c in st['srcchains'].values()))
    ws = sum(dilworth(sorted(v), leq) for v in by_src.values())
    wt = sum(dilworth(sorted(v), leq) for v in by_tgt.values())
    pk_s = h * sum(len(c) for c in st['srcchains'].values())
    pk_t = h * sum(len(c) for c in st['tgtchains'].values())
    rects = h * len(st['rects'])
    v = len(by_src)
    return dict(h=h, v=v, rects=rects, perkernel_src=pk_s, perkernel_tgt=pk_t, cross_src=ws, cross_tgt=wt,
                side_old=rects + pk_s + pk_t - 2 * v, side_new=rects + ws + wt - 2 * v)


if __name__ == "__main__":
    for h in map(int, sys.argv[1:]):
        print(run(h), flush=True)
