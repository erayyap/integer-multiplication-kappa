"""Copy-wire side cover for one kernel point q (problem on n = h-1 other points).
Elements: pairs P (source copies) and Q (target copies).  Cover every ordered (P,Q) with
P cap Q = {} exactly once by rectangles node_S x node_T with vertex-disjoint supports.
Each cell is assigned a rectangle by a deterministic rule; we verify that every rectangle
gets exactly |G||H| cells (=> exact partition).  Cost = #rectangles + chain covers."""
import sys, itertools
from collections import defaultdict
import networkx as nx

def build_tree(lo, hi):
    # nodes are (lo,hi) ranges; children by halving
    kids = {}
    def rec(lo, hi):
        if hi - lo == 1: return
        mid = (lo + hi) // 2
        kids[(lo, hi)] = ((lo, mid), (mid, hi)); rec(lo, mid); rec(mid, hi)
    rec(lo, hi); return kids

def run(n, verify=True, full=False):
    kids = build_tree(0, n)
    root = (0, n)
    def child_containing(c, x):
        a, b = kids[c]; return a if a[0] <= x < a[1] else b
    def lca(pts):
        c = root
        while c in kids:
            a, b = kids[c]
            if all(a[0] <= p < a[1] for p in pts): c = a
            elif all(b[0] <= p < b[1] for p in pts): c = b
            else: return c
        return c
    def inn(c, x): return c[0] <= x < c[1]
    def W(c): return ('W', c)
    def G(a, b): return ('G',) + tuple(sorted((a, b)))
    def rect(P, Q):
        x, y = P; z, w = Q
        c = lca((x, y, z, w)); c1, c2 = kids[c]
        side = lambda p: 0 if inn(c1, p) else 1
        sP = {side(x), side(y)}; sQ = {side(z), side(w)}
        if len(sP) == 1 and len(sQ) == 1:
            return (W(kids[c][side(x)]), W(kids[c][side(z)]))
        if len(sP) == 1:  # P inside one child, Q crossing
            s = side(x); cs, co = kids[c][s], kids[c][1 - s]
            zz = z if side(z) == s else w
            d = lca((x, y, zz)); d1, d2 = kids[d]
            if inn(d1, x) == inn(d1, y):
                dp = d1 if inn(d1, x) else d2; dz = d2 if dp == d1 else d1
                return (W(dp), G(dz, co))
            xx = x if inn(child_containing(d, zz), x) else y  # P's point on z's side
            yy = y if xx == x else x
            e = lca((xx, zz))
            return (G(child_containing(e, xx), child_containing(d, yy)), G(child_containing(e, zz), co))
        if len(sQ) == 1:  # Q inside one child, P crossing (mirror)
            s = side(z); co = kids[c][1 - s]
            xx = x if side(x) == s else y
            d = lca((z, w, xx)); d1, d2 = kids[d]
            if inn(d1, z) == inn(d1, w):
                dq = d1 if inn(d1, z) else d2; dx = d2 if dq == d1 else d1
                return (G(dx, co), W(dq))
            zz = z if inn(child_containing(d, xx), z) else w
            ww = w if zz == z else z
            e = lca((xx, zz))
            return (G(child_containing(e, xx), co), G(child_containing(e, zz), child_containing(d, ww)))
        # both crossing
        xa, xb = (x, y) if side(x) == 0 else (y, x)
        za, zb = (z, w) if side(z) == 0 else (w, z)
        e = lca((xa, za)); f = lca((xb, zb))
        return (G(child_containing(e, xa), child_containing(f, xb)),
                G(child_containing(e, za), child_containing(f, zb)))
    pairs = list(itertools.combinations(range(n), 2))
    cells = defaultdict(int); src_nodes = defaultdict(set); tgt_nodes = defaultdict(set)
    for P in pairs:
        for Q in pairs:
            if set(P) & set(Q): continue
            r = rect(P, Q); cells[r] += 1
            src_nodes[P].add(r[0]); tgt_nodes[Q].add(r[1])
    def members(node):
        if node[0] == 'W':
            c = node[1]; return [p for p in pairs if inn(c, p[0]) and inn(c, p[1])]
        a, b = node[1], node[2]
        return [p for p in pairs if (inn(a, p[0]) and inn(b, p[1])) or (inn(b, p[0]) and inn(a, p[1]))]
    def support(node):
        cs = [node[1]] if node[0] == 'W' else [node[1], node[2]]
        return set(i for c in cs for i in range(c[0], c[1]))
    if verify:
        msz = {}
        for (g, hh), k in cells.items():
            for nd in (g, hh):
                if nd not in msz: msz[nd] = len(members(nd))
            assert msz[g] * msz[hh] == k, (g, hh, k)
            assert not (support(g) & support(hh))
    def chains(nodes):
        nodes = list(nodes); S = [support(nd) for nd in nodes]
        sets = [frozenset(members(nd)) for nd in nodes]
        Gb = nx.DiGraph(); L = len(nodes)
        Gb.add_nodes_from([('l', i) for i in range(L)] + [('r', i) for i in range(L)])
        for i in range(L):
            for j in range(L):
                if i != j and sets[i] < sets[j]: Gb.add_edge(('l', i), ('r', j))
        M = nx.bipartite.hopcroft_karp_matching(Gb.to_undirected(), top_nodes=[('l', i) for i in range(L)])
        return L - len(M) // 2
    if full:
        def decompose(nodes):
            nodes = sorted(nodes, key=lambda nd: len(members(nd))); sets = [frozenset(members(nd)) for nd in nodes]
            L = len(nodes); Gb = nx.Graph(); Gb.add_nodes_from([('l', i) for i in range(L)] + [('r', i) for i in range(L)])
            for i in range(L):
                for j in range(L):
                    if sets[i] < sets[j]: Gb.add_edge(('l', i), ('r', j))
            M = nx.bipartite.hopcroft_karp_matching(Gb, top_nodes=[('l', i) for i in range(L)])
            nxt = {i: M[('l', i)][1] for i in range(L) if ('l', i) in M}
            starts = set(range(L)) - set(nxt.values()); chs = []
            for s0 in sorted(starts):
                c = [s0]
                while c[-1] in nxt: c.append(nxt[c[-1]])
                chs.append([nodes[i] for i in c])
            return chs
        mem = {}
        for nd in set(r for k in cells for r in k): mem[nd] = members(nd)
        return dict(rects=list(cells), members=mem,
                    srcchains={P: decompose(v) for P, v in src_nodes.items()},
                    tgtchains={Q: decompose(v) for Q, v in tgt_nodes.items()})
    memo = {}
    def ch(nodes):
        key = frozenset(nodes)
        if key not in memo: memo[key] = chains(nodes)
        return memo[key]
    cs = sum(ch(v) for v in src_nodes.values()); ct = sum(ch(v) for v in tgt_nodes.values())
    return dict(n=n, rects=len(cells), cells=sum(cells.values()), src_chains=cs, tgt_chains=ct,
                total=len(cells) + cs + ct)

if __name__ == "__main__":
    for n in map(int, sys.argv[1:]):
        print(run(n), flush=True)



def structure(n):
    return run(n, verify=True, full=True)
