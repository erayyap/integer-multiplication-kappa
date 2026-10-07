"""Relaxed labels: node G gets label U(A_G)=span{t_{q+P}: P subset A_G}, A_G=supp G (gather);
node H gets U(comp supp H) (scatter).  Chains need only nested supports."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys, itertools
sys.path.insert(0, f'{ROOT}/scripts')
import networkx as nx
from collections import defaultdict
import copy_sorted as cs

def width_sets(sets):
    sets = list(set(sets)); L = len(sets)   # equal supports collapse (comparable)
    Gb = nx.Graph(); Gb.add_nodes_from([('l', i) for i in range(L)] + [('r', i) for i in range(L)])
    for i in range(L):
        for j in range(L):
            if sets[i] < sets[j]: Gb.add_edge(('l', i), ('r', j))
    M = nx.bipartite.hopcroft_karp_matching(Gb, top_nodes=[('l', i) for i in range(L)])
    return L - len(M) // 2

def run(n, mod=cs):
    pairs, members = mod.members_fn(n)
    supp = {}
    def sp(nd):
        if nd not in supp: supp[nd] = frozenset(i for p in members(nd) for i in p)
        return supp[nd]
    rects = set(); src = defaultdict(set); tgt = defaultdict(set)
    for P in pairs:
        for Q in pairs:
            if set(P) & set(Q): continue
            r = mod.rule(P, Q); rects.add(r); src[P].add(sp(r[0])); tgt[Q].add(sp(r[1]))
    a = sum(width_sets(v) for v in src.values()); b = sum(width_sets(v) for v in tgt.values())
    return dict(n=n, rects=len(rects), src_chains=a, tgt_chains=b, total=len(rects) + a + b)

if __name__ == "__main__":
    for n in map(int, sys.argv[1:]): print(run(n), flush=True)
