"""Sorted-order copy-wire side cover for one kernel point (n other points 0..n-1).
Cells (P,Q), P=(x1<x2), Q=(z1<z2), disjoint; rectangle chosen by the order pattern."""
import sys, itertools
from collections import defaultdict
import networkx as nx

def members_fn(n):
    pairs = list(itertools.combinations(range(n), 2))
    def members(nd):
        t = nd[0]
        if t == 'Wpre': z = nd[1]; return [p for p in pairs if p[1] < z]
        if t == 'rowfrom': x = nd[1]; return [p for p in pairs if p[0] == x]
        if t == 'rowsuf': x, z = nd[1:]; return [p for p in pairs if p[0] == x and p[1] > z]
        if t == 'colpre': z, x = nd[1:]; return [p for p in pairs if p[1] == x and p[0] < z]
        if t == 'colmid': z, x = nd[1:]; return [p for p in pairs if p[1] == x and z < p[0]]
        raise ValueError(nd)
    return pairs, members

def rule(P, Q):
    x1, x2 = P; z1, z2 = Q
    if x2 < z1: return (('Wpre', z1), ('rowfrom', z1))
    if z2 < x1: return (('rowfrom', x1), ('Wpre', x1))
    if x1 < z1 < x2 < z2: return (('colpre', z1, x2), ('rowsuf', z1, x2))
    if z1 < x1 < z2 < x2: return (('rowsuf', x1, z2), ('colpre', x1, z2))
    if x1 < z1 < z2 < x2: return (('rowsuf', x1, z2), ('colmid', x1, z2))
    if z1 < x1 < x2 < z2: return (('colmid', z1, x2), ('rowsuf', z1, x2))
    raise ValueError

def width(nodes, members):
    nodes = list(nodes); sets = [frozenset(members(nd)) for nd in nodes]; L = len(nodes)
    Gb = nx.Graph(); Gb.add_nodes_from([('l', i) for i in range(L)] + [('r', i) for i in range(L)])
    for i in range(L):
        for j in range(L):
            if sets[i] < sets[j]: Gb.add_edge(('l', i), ('r', j))
    M = nx.bipartite.hopcroft_karp_matching(Gb, top_nodes=[('l', i) for i in range(L)])
    return L - len(M) // 2

def run(n, rule=rule, verify=True):
    pairs, members = members_fn(n)
    cells = defaultdict(int); src = defaultdict(set); tgt = defaultdict(set)
    for P in pairs:
        for Q in pairs:
            if set(P) & set(Q): continue
            r = rule(P, Q); cells[r] += 1; src[P].add(r[0]); tgt[Q].add(r[1])
    if verify:
        for (g, hh), k in cells.items():
            mg, mh = members(g), members(hh)
            assert len(mg) * len(mh) == k, (g, hh, k)
            sg = set(i for p in mg for i in p); sh = set(i for p in mh for i in p)
            assert not (sg & sh)
    cs = sum(width(v, members) for v in src.values()); ct = sum(width(v, members) for v in tgt.values())
    return dict(n=n, rects=len(cells), cells=sum(cells.values()), src_chains=cs, tgt_chains=ct, total=len(cells) + cs + ct)

if __name__ == "__main__":
    for n in map(int, sys.argv[1:]): print(run(n), flush=True)


def decompose(nodes, members):
    nodes = sorted(nodes, key=lambda nd: len(members(nd))); sets = [frozenset(members(nd)) for nd in nodes]
    L = len(nodes); Gb = nx.Graph(); Gb.add_nodes_from([('l', i) for i in range(L)] + [('r', i) for i in range(L)])
    for i in range(L):
        for j in range(L):
            if sets[i] < sets[j]: Gb.add_edge(('l', i), ('r', j))
    M = nx.bipartite.hopcroft_karp_matching(Gb, top_nodes=[('l', i) for i in range(L)])
    nxt = {i: M[('l', i)][1] for i in range(L) if ('l', i) in M}
    chs = []
    for s0 in sorted(set(range(L)) - set(nxt.values())):
        c = [s0]
        while c[-1] in nxt: c.append(nxt[c[-1]])
        chs.append([nodes[i] for i in c])
    return chs


def structure(n, rule=rule):
    pairs, members = members_fn(n)
    cells = defaultdict(int); src = defaultdict(set); tgt = defaultdict(set)
    for P in pairs:
        for Q in pairs:
            if set(P) & set(Q): continue
            r = rule(P, Q); cells[r] += 1; src[P].add(r[0]); tgt[Q].add(r[1])
    mem = {nd: members(nd) for r in cells for nd in r}
    for (g, hh), k in cells.items(): assert len(mem[g]) * len(mem[hh]) == k
    return dict(rects=list(cells), members=mem,
                srcchains={P: decompose(v, members) for P, v in src.items()},
                tgtchains={Q: decompose(v, members) for Q, v in tgt.items()})
