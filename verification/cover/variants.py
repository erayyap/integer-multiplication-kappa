import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys, itertools
from collections import defaultdict
sys.path.insert(0, f'{ROOT}/scripts')
from copy_sorted import width
def mk(optA, optB):
    def rule(P, Q):
        x1, x2 = P; z1, z2 = Q
        if x2 < z1: return ('sepG', z1)
        if z2 < x1: return ('sepH', x1)
        if x1 < z1 < x2 < z2: return ('PQPQ', z1, x2)
        if z1 < x1 < z2 < x2: return ('QPQP', x1, z2)
        if x1 < z1 < z2 < x2: return ('PQQP-A', x1, z2) if optA else ('PQQP-B', x2, z1)
        if z1 < x1 < x2 < z2: return ('QPPQ-A', z1, x2) if optB else ('QPPQ-B', x1, z2)
    return rule
def evaluate(n, rule):
    pairs = list(itertools.combinations(range(n), 2))
    G = defaultdict(set); H = defaultdict(set); cnt = defaultdict(int)
    for P in pairs:
        for Q in pairs:
            if set(P) & set(Q): continue
            k = rule(P, Q); G[k].add(P); H[k].add(Q); cnt[k] += 1
    for k in cnt: assert len(G[k]) * len(H[k]) == cnt[k], k
    src = defaultdict(set); tgt = defaultdict(set)
    for k in cnt:
        g = frozenset(G[k]); h = frozenset(H[k])
        for P in g: src[P].add(g)
        for Q in h: tgt[Q].add(h)
    ident = lambda s: s
    cs = sum(width(v, ident) for v in src.values()); ct = sum(width(v, ident) for v in tgt.values())
    return len(cnt), cs, ct, len(cnt) + cs + ct
for n in map(int, sys.argv[1:]):
    for a in (True, False):
        for b in (True, False):
            print(n, 'PQQP', 'A' if a else 'B', 'QPPQ', 'A' if b else 'B', evaluate(n, mk(a, b)), flush=True)
