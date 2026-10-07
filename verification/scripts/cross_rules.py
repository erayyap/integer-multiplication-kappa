"""cross-kernel widths (rot order) for rule variants given as cell -> key."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys, itertools
sys.path.insert(0, f'{ROOT}/scripts'); sys.path.insert(0, f'{ROOT}/cover')
import cross_width as cw
from collections import defaultdict
def structure(n, rule):
    pairs = list(itertools.combinations(range(n), 2)); G = defaultdict(set); H = defaultdict(set); cnt = defaultdict(int)
    for P in pairs:
        for Q in pairs:
            if set(P) & set(Q): continue
            k = rule(P, Q); G[k].add(P); H[k].add(Q); cnt[k] += 1
    for k in cnt: assert len(G[k]) * len(H[k]) == cnt[k], k
    return G, H
def run(h, rule, order):
    n = h - 1; G, H = structure(n, rule)
    ind = lambda t: [1 if i in t else 0 for i in range(h)]
    TN = {}
    for q in range(h):
        o = order(q, h)
        for k in G: TN[(q, 'g', k)] = frozenset(tuple(sorted((q, o[a], o[b]))) for a, b in G[k])
        for k in H: TN[(q, 'h', k)] = frozenset(tuple(sorted((q, o[a], o[b]))) for a, b in H[k])
    B = {k: cw.rref([ind(t) for t in v], h) for k, v in TN.items()}
    c = {}
    def leq(a, b):
        if (a, b) not in c: c[(a, b)] = TN[a] <= TN[b] or all(cw.inspan(ind(t), B[b]) for t in TN[a])
        return c[(a, b)]
    key = lambda k: (len(B[k]), str(k))
    tot = 0
    for role in ('g', 'h'):
        by = defaultdict(set)
        for k, v in TN.items():
            if k[1] == role:
                for t in v: by[t].add(k)
        tot += sum(cw.dilworth(sorted(v, key=key), leq) for v in by.values())
    return h * len(G), tot
if __name__ == "__main__":
    from variants import mk
    rot = lambda q, h: [(q + 1 + i) % h for i in range(h - 1)]
    h = int(sys.argv[1])
    for a in (True, False):
        for b in (True, False):
            r, ch = run(h, mk(a, b), rot); print(h, a, b, 'rects', r, 'chains', ch, 'sum', r + ch, flush=True)
