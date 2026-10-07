"""rot order, rule variant chosen per kernel by a function of q."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys, io, contextlib
sys.path.insert(0, f'{ROOT}/scripts'); sys.path.insert(0, f'{ROOT}/cover')
import cross_width as cw, cross_rules as cr
from collections import defaultdict
with contextlib.redirect_stdout(io.StringIO()):
    from variants import mk
def run(h, choose):
    n = h - 1; S = {}
    for v in ('AA', 'AB', 'BA', 'BB'): S[v] = cr.structure(n, mk(v[0] == 'A', v[1] == 'A'))
    ind = lambda t: [1 if i in t else 0 for i in range(h)]
    TN = {}; R = 0
    for q in range(h):
        o = [(q + 1 + i) % h for i in range(h - 1)]; v = choose(q, h); G, H = S[v]; R += len(G)
        for k in G: TN[(q, 'g', k)] = frozenset(tuple(sorted((q, o[a], o[b]))) for a, b in G[k])
        for k in H: TN[(q, 'h', k)] = frozenset(tuple(sorted((q, o[a], o[b]))) for a, b in H[k])
    B = {k: cw.rref([ind(t) for t in v], h) for k, v in TN.items()}
    c = {}
    def leq(a, b):
        if (a, b) not in c: c[(a, b)] = TN[a] <= TN[b] or all(cw.inspan(ind(t), B[b]) for t in TN[a])
        return c[(a, b)]
    tot = R
    for role in ('g', 'h'):
        by = defaultdict(set)
        for k, v in TN.items():
            if k[1] == role:
                for t in v: by[t].add(k)
        tot += sum(cw.dilworth(sorted(v, key=lambda k: (len(B[k]), str(k))), leq) for v in by.values())
    return tot
if __name__ == "__main__":
    h = int(sys.argv[1])
    C = {'allAA': lambda q, h: 'AA', 'parAA_AB': lambda q, h: 'AA' if q % 2 == 0 else 'AB',
         'parAA_BB': lambda q, h: 'AA' if q % 2 == 0 else 'BB', 'parAB_BA': lambda q, h: 'AB' if q % 2 == 0 else 'BA',
         'halfAA_BB': lambda q, h: 'AA' if q < h // 2 else 'BB', 'mod3': lambda q, h: ('AA', 'AB', 'BB')[q % 3]}
    for k, f in C.items(): print(h, k, run(h, f), flush=True)
