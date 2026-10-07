"""Greedy merging of side wires (rects) that share a source triple set (or a target triple set).

Rect (G, H) and (G, H') with the same triple set G -> one wire (G, H u H'), a new target node.
Valid: G x (H u H') covers the same cells; G _|_ H u H'.  Cost = #rects + sum of cross-kernel
Dilworth widths (sources and targets), nodes identified by triple sets (equal sets = one gate).
Only accepts a merge if the total strictly decreases.  Also checks nondegeneracy of merged spans.
"""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys, itertools, random
from collections import defaultdict
sys.path.insert(0, f'{ROOT}/scripts')
import cross_width as cw, copy_sorted as cs

h = int(sys.argv[1]); n = h - 1; seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
p = cw.p
inv9 = pow(9, p - 2, p)
import os
VAR = os.environ.get('VAR', 'AA')
if VAR == 'AA':
    st = cs.structure(n); mem = st['members']; RECTS = st['rects']
else:
    sys.path.insert(0, f'{ROOT}/cover')
    import io, contextlib, cross_rules
    with contextlib.redirect_stdout(io.StringIO()):
        from variants import mk
    Gd, Hd = cross_rules.structure(n, mk(VAR[0] == 'A', VAR[1] == 'A'))
    mem = {}; RECTS = []
    for k in Gd:
        mem[('g', k)] = sorted(Gd[k]); mem[('h', k)] = sorted(Hd[k]); RECTS.append((('g', k), ('h', k)))
rot = lambda q: [(q + 1 + i) % h for i in range(h - 1)]
ind = lambda t: [1 if i in t else 0 for i in range(h)]
rects = []
for q in range(h):
    o = rot(q)
    T = lambda nd: frozenset(tuple(sorted((q, o[a], o[b]))) for a, b in mem[nd])
    for g, hh in RECTS: rects.append([T(g), T(hh)])
B = {}
def basis(s):
    if s not in B: B[s] = cw.rref([ind(t) for t in s], h)
    return B[s]
C = {}
def leq(a, b):
    if (a, b) not in C: C[(a, b)] = a <= b or all(cw.inspan(ind(t), basis(b)) for t in a)
    return C[(a, b)]
def width(nodes):
    items = sorted(nodes, key=lambda s: (len(basis(s)), sorted(s)))
    return cw.dilworth(items, leq)
def nondeg(s):
    bs = [b for _, b in basis(s)]
    G = [[(sum(x * y for x, y in zip(u, w)) - inv9 * sum(u) * sum(w)) % p for w in bs] for u in bs]
    return len(cw.rref(G, len(bs))) == len(bs)
def nodes_of(rs, role):
    by = defaultdict(set)
    for r in rs:
        for t in r[role]: by[t].add(r[role])
    return by
def total(rs):
    s = len(rs)
    for role in (0, 1):
        s += sum(width(v) for v in nodes_of(rs, role).values())
    return s
base = total(rects); print('h', h, 'rects', len(rects), 'total', base, flush=True)
random.seed(seed)
improved = True; cur = base
while improved:
    improved = False
    for role in (0, 1):
        other = 1 - role
        groups = defaultdict(list)
        for i, r in enumerate(rects):
            if r is not None: groups[r[role]].append(i)
        cands = [(i, j) for g in groups.values() for i, j in itertools.combinations(g, 2)]
        random.shuffle(cands)
        for i, j in cands:
            if rects[i] is None or rects[j] is None or rects[i][role] != rects[j][role]: continue
            U = rects[i][other] | rects[j][other]
            if not nondeg(U): continue
            # local delta: elements touching the other-role sets
            aff = set(rects[i][other]) | set(rects[j][other])
            live = [r for r in rects if r is not None]
            by = nodes_of(live, other)
            before = sum(width(by[t]) for t in aff)
            new = [r for k, r in enumerate(rects) if r is not None and k not in (i, j)]
            m = [None, None]; m[role] = rects[i][role]; m[other] = U; new.append(m)
            by2 = nodes_of(new, other)
            after = sum(width(by2[t]) for t in aff)
            delta = -1 + after - before
            if delta < 0:
                rects[i] = m; rects[j] = None; cur += delta; improved = True
                print('merge role', role, 'delta', delta, 'total', cur, flush=True)
    rects = [r for r in rects if r is not None]
print('final total', cur, 'check', total(rects), 'rects', len(rects), flush=True)
