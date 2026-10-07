"""Explicit check of one stage-1 invocation of the bit network with copy-wire side covers.
Stage 1: A = ground field, B = 0, so every gate label is B(x)F + P(x)U = U (U subspace of F=Q^h)
or B(x)<t> = 0.  Form on F: I - J/9, t_S = indicator of S.
Builds every wire (data X,Y; centers; side wires s_{G,H}; source copies; target accumulators),
an explicit time-ordered gate list with labels, and checks
 (1) F2: for random inputs and random initial scratch values, y_T -> y_T + x_T, x and all
     scratch wires restored;
 (2) every wire's label sequence (source label, gates in time order, sink label) is nested
     and the only decreasing edges are centre 3->4 (loss h each, as in the paper);
 (3) every gate label is nondegenerate.
Then prints the wire count."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys, itertools, random
from fractions import Fraction as Fr
from collections import defaultdict
import numpy as np
import networkx as nx
sys.path.insert(0, f'{ROOT}/scripts')
import importlib; cc = importlib.import_module(sys.argv[2] if len(sys.argv) > 2 else "copy_cover")

h = int(sys.argv[1]); n = h - 1
p = 1000003
inv9 = pow(9, p - 2, p)
Gm = (np.eye(h, dtype=np.int64) - inv9 * np.ones((h, h), dtype=np.int64)) % p

def rank(M):
    M = np.array(M, dtype=np.int64) % p
    if M.size == 0: return 0
    r = 0; rows, cols = M.shape
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i, c]), None)
        if piv is None: continue
        M[[r, piv]] = M[[piv, r]]; M[r] = M[r] * pow(int(M[r, c]), p - 2, p) % p
        for i in range(rows):
            if i != r and M[i, c]: M[i] = (M[i] - M[i, c] * M[r]) % p
        r += 1
        if r == rows: break
    return r
def nullspace(M):
    M = np.array(M, dtype=np.int64) % p; rows, cols = M.shape; r = 0; piv = []
    for c in range(cols):
        pv = next((i for i in range(r, rows) if M[i, c]), None)
        if pv is None: continue
        M[[r, pv]] = M[[pv, r]]; M[r] = M[r] * pow(int(M[r, c]), p - 2, p) % p
        for i in range(rows):
            if i != r and M[i, c]: M[i] = (M[i] - M[i, c] * M[r]) % p
        piv.append(c); r += 1
        if r == rows: break
    out = []
    for f in [c for c in range(cols) if c not in piv]:
        x = np.zeros(cols, dtype=np.int64); x[f] = 1
        for i, c in enumerate(piv): x[c] = (-M[i, f]) % p
        out.append(x)
    return np.array(out, dtype=np.int64).reshape(len(out), cols)
def perp(U):
    if len(U) == 0: return np.eye(h, dtype=np.int64)
    return nullspace(U @ Gm % p)
def contained(U, V):
    if len(U) == 0: return True
    if len(V) == 0: return rank(U) == 0
    return rank(np.vstack([V, U])) == rank(V)
def nondeg(U):
    return len(U) == 0 or rank(U @ Gm @ U.T % p) == rank(U)

ZERO = np.zeros((0, h), dtype=np.int64); FULL = np.eye(h, dtype=np.int64)
T = list(itertools.combinations(range(h), 3)); v = len(T)
ind = {t: np.array([1 if i in t else 0 for i in range(h)], dtype=np.int64) for t in T}

# per kernel point q: copy_cover structure on the other points
def family(q):
    others = [i for i in range(h) if i != q]
    st = cc.structure(n)   # rects (srcnode, tgtnode), node members as index pairs
    def trip(P): return tuple(sorted((q, others[P[0]], others[P[1]])))
    return st, trip

wires = {}   # name -> list of (time, label, kind)
def touch(w, t, U): wires.setdefault(w, []).append((t, U))
gates = []   # (time, label, list of (target_wire, [source wires]))  F2: target ^= xor(sources)

def add_gate(t, U, ops):
    gates.append((t, U, ops))
    for tgt, srcs in ops:
        touch(tgt, t, U)
        for s in srcs: touch(s, t, U)

def span(trips): return np.array([ind[x] for x in trips], dtype=np.int64)

# time layout (fractional times inside the paper's rows)
# 0: y -= b (BF);  0.5: b -= s at H (BF); 1: centre subtract scatter (BF); 1.5: garbage s -= a at G (BF)
# 2: own copies a += x (span t_X); 2.x: source nodes bottom-up s += a (span G)
# 3: centre gather (AF); 4: centre scatter (BF); 5.x: target nodes top-down b += s (V_H^perp)
# 5.9: own' y += b (t_T^perp); 6: centre undo gather (AF); 7.x: restores (AF)
import os
# ---- cross-kernel chains (check_cross): copy/accumulator chains may visit nodes of all 3 kernels
import cross_width as cw
st = cc.structure(n); rects0, members = st['rects'], st['members']
srcn = sorted({g for g, _ in rects0}, key=str); tgtn = sorted({hh for _, hh in rects0}, key=str)
ORDER = os.environ.get('ORDER', 'global')
def tripq(q, P):
    others = [i for i in range(h) if i != q] if ORDER == 'global' else [(q + 1 + i) % h for i in range(h - 1)]
    return tuple(sorted((q, others[P[0]], others[P[1]])))
NT = {}
for q in range(h):
    for nd in set(srcn) | set(tgtn): NT[(q, nd)] = [tripq(q, P) for P in members[nd]]
spanrank = {k: rank(span(vv)) for k, vv in NT.items()}
basis = {k: cw.rref([[int(x) for x in ind[t]] for t in vv], h) for k, vv in NT.items()}
gidx = {k: i for i, k in enumerate(sorted(NT, key=lambda k: (spanrank[k], str(k))))}
_c = {}
def leq(a, b):
    if (a, b) not in _c:
        _c[(a, b)] = set(NT[a]) <= set(NT[b]) or all(cw.inspan([int(x) for x in ind[t]], basis[b]) for t in NT[a])
    return _c[(a, b)]
def chains(items):
    items = sorted(items, key=lambda k: gidx[k]); L = len(items)
    Gb = nx.Graph(); Gb.add_nodes_from([('l', i) for i in range(L)] + [('r', i) for i in range(L)])
    for i in range(L):
        for j in range(i + 1, L):
            if leq(items[i], items[j]): Gb.add_edge(('l', i), ('r', j))
    M = nx.bipartite.hopcroft_karp_matching(Gb, top_nodes=[('l', i) for i in range(L)])
    nxt = {i: M[('l', i)][1] for i in range(L) if ('l', i) in M}
    out = []
    for s0 in sorted(set(range(L)) - set(nxt.values())):
        c = [s0]
        while c[-1] in nxt: c.append(nxt[c[-1]])
        out.append([items[i] for i in c])
    return out
by_src = defaultdict(list); by_tgt = defaultdict(list)
for q in range(h):
    for nd in srcn:
        for t in NT[(q, nd)]: by_src[t].append((q, nd))
    for nd in tgtn:
        for t in NT[(q, nd)]: by_tgt[t].append((q, nd))
src_of = {}; tgt_of = {}; acopies = []; bcopies = []
for t, nds in by_src.items():
    for k, ch in enumerate(chains(nds)):
        w = ('X', t) if k == 0 else ('a', t, k)
        if k: acopies.append((t, w))
        for nd in ch: src_of[(nd, t)] = w
for t, nds in by_tgt.items():
    for k, ch in enumerate(chains(nds)):
        w = ('Y', t) if k == 0 else ('b', t, k)
        if k: bcopies.append((t, w))
        for nd in ch: tgt_of[(nd, t)] = w
byG = defaultdict(list); byH = defaultdict(list)
for q in range(h):
    for r in rects0: byG[(q, r[0])].append((q, r)); byH[(q, r[1])].append((q, r))
sw = {qr: ('s',) + qr for qr in [(q, r) for q in range(h) for r in rects0]}
side_count = len(sw); copy_count = len(acopies) + len(bcopies)
SRC = sorted([(q, g) for q in range(h) for g in srcn], key=lambda k: gidx[k])
TGT = sorted([(q, g) for q in range(h) for g in tgtn], key=lambda k: -gidx[k])
def glab(k): return span(NT[k])
def hlab(k): return perp(span(NT[k]))
for t, w in bcopies: add_gate(0, ZERO, [(('Y', t), [w])])
for i, H in enumerate(TGT):
    add_gate(0.5 + i * 1e-7, ZERO, [(tgt_of[(H, t)], [sw[r] for r in byH[H]]) for t in NT[H]])
for i, G in enumerate(SRC):
    add_gate(1.5 + i * 1e-7, ZERO, [(sw[r], [src_of[(G, t)] for t in NT[G] if src_of[(G, t)][0] == 'a']) for r in byG[G]])
for t, w in acopies: add_gate(2, span([t]), [(w, [('X', t)])])
for i, G in enumerate(SRC):
    add_gate(2.1 + i * 1e-7, glab(G), [(sw[r], [src_of[(G, t)] for t in NT[G]]) for r in byG[G]])
for i, H in enumerate(TGT):
    add_gate(5.1 + i * 1e-7, hlab(H), [(tgt_of[(H, t)], [sw[r] for r in byH[H]]) for t in NT[H]])
for t, w in bcopies: add_gate(5.9, perp(span([t])), [(('Y', t), [w])])
for i, H in enumerate(TGT):
    add_gate(7.1 + i * 1e-7, FULL, [(tgt_of[(H, t)], [sw[r] for r in byH[H]]) for t in NT[H] if tgt_of[(H, t)][0] == 'b'])
for i, G in enumerate(SRC):
    add_gate(7.2 + i * 1e-7, FULL, [(sw[r], [src_of[(G, t)] for t in NT[G]]) for r in byG[G]])
for t, w in acopies: add_gate(7.3, FULL, [(w, [('X', t)])])
for i, G in enumerate(SRC):
    add_gate(7.4 + i * 1e-7, FULL, [(sw[r], [src_of[(G, t)] for t in NT[G] if src_of[(G, t)][0] == 'a']) for r in byG[G]])
for i, H in enumerate(TGT):
    add_gate(7.5 + i * 1e-7, FULL, [(tgt_of[(H, t)], [sw[r] for r in byH[H]]) for t in NT[H] if tgt_of[(H, t)][0] == 'b'])
# centres
for i in range(h):
    Ci = ('C', i)
    add_gate(1, ZERO, [(('Y', t), [Ci]) for t in T if i in t])
    add_gate(3, FULL, [(Ci, [('X', t) for t in T if i in t])])
    add_gate(4, ZERO, [(('Y', t), [Ci]) for t in T if i in t])
    add_gate(6, FULL, [(Ci, [('X', t) for t in T if i in t])])
for t in T: touch(('X', t), -1, None); touch(('Y', t), -1, None)
gates.sort(key=lambda g: g[0])
import os
DUAL = os.environ.get('DUAL') == '1'
if DUAL:
    # stage 2: inverse network (gates reversed; F2 XOR gates are involutions), logical banks
    # exchanged (rename X<->Y), every gate label replaced by its orthogonal complement in F.
    def ren(w): return ('Y', w[1]) if w[0] == 'X' else (('X', w[1]) if w[0] == 'Y' else w)
    gates = [(10 - t, perp(U), [(ren(tg), [ren(x) for x in srcs]) for tg, srcs in ops]) for t, U, ops in reversed(gates)]
    wires = {}
    for t, U, ops in gates:
        for tg, srcs in ops:
            touch(tg, t, U)
            for x in srcs: touch(x, t, U)
    for t in T: touch(('X', t), -1, None); touch(('Y', t), -1, None)
# (1) F2 simulation
names = sorted(wires, key=str); idx = {w: i for i, w in enumerate(names)}
rng = np.random.default_rng(1)
for trial in range(3):
    val = rng.integers(0, 2, len(names)); init = val.copy()
    for _, _, ops in gates:
        for tgt, srcs in ops:
            if srcs: val[idx[tgt]] ^= np.bitwise_xor.reduce(val[[idx[s] for s in srcs]])
    for w in names:
        upd, ctl = ('X', 'Y') if DUAL else ('Y', 'X')
        if w[0] == upd: assert val[idx[w]] == init[idx[w]] ^ init[idx[(ctl, w[1])]], w
        else: assert val[idx[w]] == init[idx[w]], w
print(f"h={h}: F2 net map {'x += y (stage 2, dual network)' if DUAL else 'y += x'}, all scratch restored (3 random trials)")
# (2) label chains.  source/sink labels: X: <t> -> F ; Y: 0 -> t^perp ; scratch: 0 -> F
lab_ok = 0; dec = 0
for w, seq in wires.items():
    seq = sorted([s for s in seq if s[1] is not None], key=lambda s: s[0])
    if w[0] == 'X': src, snk = span([w[1]]), FULL
    elif w[0] == 'Y': src, snk = ZERO, perp(span([w[1]]))
    else: src, snk = ZERO, FULL
    labs = [src] + [U for _, U in seq] + [snk]
    # collapse identical consecutive labels quickly
    for a, b in zip(labs, labs[1:]):
        if not contained(a, b):
            assert w[0] == 'C' and rank(a) == h and rank(b) == 0, (w, rank(a), rank(b))
            dec += 1
    lab_ok += 1
assert dec == 2 * h * 0 + h, dec  # centre 3->4 only (one per centre wire)
print(f"label chains nested on all {lab_ok} wires; decreasing edges = {dec} (centre 3->4 only)")
# (3) nondegenerate gate labels
seen = set(); bad = 0
for _, U, _ in gates:
    key = U.tobytes()
    if key in seen: continue
    seen.add(key); bad += (not nondeg(U))
assert bad == 0
print(f"all {len(seen)} distinct gate labels nondegenerate")
print(f"wires: data {2*v}, centre {h}, side {side_count}, copies {copy_count}, scratch total {len(wires)-2*v}")
