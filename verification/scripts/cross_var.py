"""cross-kernel widths for the sorted cover with per-kernel point orders."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys
sys.path.insert(0, f'{ROOT}/scripts')
import cross_width as cw, copy_sorted as cs
from collections import defaultdict
def run(h, order, rule=None):
    n = h - 1
    st = cs.structure(n) if rule is None else cs.structure(n, rule)
    mem = st['members']; srcn = {g for g, _ in st['rects']}; tgtn = {x for _, x in st['rects']}
    ind = lambda t: [1 if i in t else 0 for i in range(h)]
    TN = {}
    for q in range(h):
        o = order(q, h)
        for nd in srcn | tgtn: TN[(q, nd)] = frozenset(tuple(sorted((q, o[a], o[b]))) for a, b in mem[nd])
    B = {k: cw.rref([ind(t) for t in v], h) for k, v in TN.items()}
    c = {}
    def leq(a, b):
        if (a, b) not in c: c[(a, b)] = TN[a] <= TN[b] or all(cw.inspan(ind(t), B[b]) for t in TN[a])
        return c[(a, b)]
    rk = {k: len(B[k]) for k in TN}
    key = lambda k: (rk[k], str(k))
    tot = 0
    for nodes in (srcn, tgtn):
        by = defaultdict(set)
        for q in range(h):
            for nd in nodes:
                for t in TN[(q, nd)]: by[t].add((q, nd))
        tot += sum(cw.dilworth(sorted(v, key=key), leq) for v in by.values())
    return tot
ORD = {
 'global': lambda q, h: [i for i in range(h) if i != q],
 'rot': lambda q, h: [(q + 1 + i) % h for i in range(h - 1)],
 'rotrev': lambda q, h: [(q - 1 - i) % h for i in range(h - 1)],
 'altrev': lambda q, h: [i for i in range(h) if i != q][::(1 if q % 2 == 0 else -1)],
 'reverse': lambda q, h: [i for i in range(h) if i != q][::-1],
 'aboveFirst': lambda q, h: [i for i in range(q + 1, h)] + [i for i in range(q)],
 'belowRevFirst': lambda q, h: [i for i in range(q - 1, -1, -1)] + [i for i in range(q + 1, h)],
}
if __name__ == "__main__":
    h = int(sys.argv[1])
    for name, o in ORD.items(): print(h, name, run(h, o), flush=True)
