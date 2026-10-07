"""copy_sorted cover with relaxed (support) labels: A_G = supp G, B_H = supp H; chains = Dilworth
decomposition of the support posets (equal supports comparable).  Test module for check_copy."""
import networkx as nx
import copy_sorted as cs

def chains_by(nodes, key, increasing=True):
    nodes = sorted(nodes, key=lambda nd: (len(key[nd]), str(nd)))
    if not increasing: nodes = nodes  # chains listed small->large support; caller orders time
    L = len(nodes); sets = [key[nd] for nd in nodes]
    Gb = nx.Graph(); Gb.add_nodes_from([('l', i) for i in range(L)] + [('r', i) for i in range(L)])
    for i in range(L):
        for j in range(i + 1, L):
            if sets[i] <= sets[j]: Gb.add_edge(('l', i), ('r', j))
    M = nx.bipartite.hopcroft_karp_matching(Gb, top_nodes=[('l', i) for i in range(L)])
    nxt = {i: M[('l', i)][1] for i in range(L) if ('l', i) in M}
    out = []
    for s0 in sorted(set(range(L)) - set(nxt.values())):
        c = [s0]
        while c[-1] in nxt: c.append(nxt[c[-1]])
        out.append([nodes[i] for i in c])
    return out

def structure(n):
    st = cs.structure(n); mem = st['members']
    supp = {nd: frozenset(i for p in m for i in p) for nd, m in mem.items()}
    src = {}; tgt = {}
    for g, hh in st['rects']:
        for P in mem[g]: src.setdefault(P, set()).add(g)
        for Q in mem[hh]: tgt.setdefault(Q, set()).add(hh)
    st['gsupp'] = {g: supp[g] for g, _ in st['rects']}
    st['hsupp'] = {hh: supp[hh] for _, hh in st['rects']}
    st['srcchains'] = {P: chains_by(v, supp) for P, v in src.items()}
    st['tgtchains'] = {Q: chains_by(v, supp) for Q, v in tgt.items()}
    return st
