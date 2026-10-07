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
side_count = 0; copy_count = 0
for q in range(h):
    st, trip = family(q)
    rects, members = st['rects'], st['members']
    srcchains, tgtchains = st['srcchains'], st['tgtchains']   # element -> list of chains (lists of nodes, small->large)
    # wires for chains: first source chain of X uses the data wire if q == min(X); likewise targets
    def swire(P, k):
        X = trip(P)
        return ('X', X) if (k == 0 and q == min(X)) else ('a', q, P, k)
    def twire(Q, k):
        Y = trip(Q)
        return ('Y', Y) if (k == 0 and q == min(Y)) else ('b', q, Q, k)
    src_of = {}  # (P,node) -> wire
    for P, chs in srcchains.items():
        for k, ch in enumerate(chs):
            for nd in ch: src_of[(P, nd)] = swire(P, k)
    tgt_of = {}
    for Q, chs in tgtchains.items():
        for k, ch in enumerate(chs):
            for nd in ch: tgt_of[(Q, nd)] = twire(Q, k)
    gs, hs = st.get('gsupp'), st.get('hsupp')   # optional relaxed labels: point sets A_G, B_H
    def UA(A): return span([trip(P) for P in itertools.combinations(sorted(A), 2)])
    def sp(nd): return span([trip(P) for P in members[nd]])
    def glab(nd): return UA(gs[nd]) if gs else sp(nd)
    def hlab(nd): return perp(UA(hs[nd])) if hs else perp(sp(nd))
    srcnodes = sorted({g for g, _ in rects}, key=lambda nd: len(gs[nd]) if gs else len(members[nd]))
    tgtnodes = sorted({hh for _, hh in rects}, key=lambda nd: -len(hs[nd]) if hs else -len(members[nd]))
    if gs:
        for g in srcnodes: assert all(set(P) <= gs[g] for P in members[g])
    if hs:
        for hh in tgtnodes: assert all(set(Q) <= hs[hh] for Q in members[hh])
    byG = defaultdict(list); byH = defaultdict(list)
    for r in rects: byG[r[0]].append(r); byH[r[1]].append(r)
    sw = {r: ('s', q, r) for r in rects}; side_count += len(rects)
    acopies = [(P, k) for P, chs in srcchains.items() for k in range(len(chs)) if swire(P, k)[0] == 'a']
    bcopies = [(Q, k) for Q, chs in tgtchains.items() for k in range(len(chs)) if twire(Q, k)[0] == 'b']
    copy_count += len(acopies) + len(bcopies)
    # 0: y -= b
    for Q, k in bcopies: add_gate(0, ZERO, [(('Y', trip(Q)), [twire(Q, k)])])
    # 0.5: b -= s at H, label BF (=0 at stage 1)
    for i, H in enumerate(tgtnodes):
        add_gate(0.5 + i * 1e-6, ZERO, [(tgt_of[(Q, H)], [sw[r] for r in byH[H]]) for Q in members[H]])
    # 1.5 garbage: s -= a for copy members
    for i, G in enumerate(srcnodes):
        ops = [(sw[r], [src_of[(P, G)] for P in members[G] if src_of[(P, G)][0] == 'a']) for r in byG[G]]
        add_gate(1.5 + i * 1e-6, ZERO, ops)
    # 2: own copies
    for P, k in acopies: add_gate(2, span([trip(P)]), [(swire(P, k), [('X', trip(P))])])
    # 2.x: source nodes
    for i, G in enumerate(srcnodes):
        add_gate(2.1 + i * 1e-6, glab(G), [(sw[r], [src_of[(P, G)] for P in members[G]]) for r in byG[G]])
    # 5.x: target nodes, decreasing span
    for i, H in enumerate(tgtnodes):
        add_gate(5.1 + i * 1e-6, hlab(H), [(tgt_of[(Q, H)], [sw[r] for r in byH[H]]) for Q in members[H]])
    for Q, k in bcopies: add_gate(5.9, perp(span([trip(Q)])), [(('Y', trip(Q)), [twire(Q, k)])])
    # 7.x restores (all AF)
    for i, H in enumerate(tgtnodes):
        add_gate(7.1 + i * 1e-6, FULL, [(tgt_of[(Q, H)], [sw[r] for r in byH[H]]) for Q in members[H] if tgt_of[(Q, H)][0] == 'b'])
    for i, G in enumerate(srcnodes):
        add_gate(7.2 + i * 1e-6, FULL, [(sw[r], [src_of[(P, G)] for P in members[G]]) for r in byG[G]])
    for P, k in acopies: add_gate(7.3, FULL, [(swire(P, k), [('X', trip(P))])])
    for i, G in enumerate(srcnodes):
        add_gate(7.4 + i * 1e-6, FULL, [(sw[r], [src_of[(P, G)] for P in members[G] if src_of[(P, G)][0] == 'a']) for r in byG[G]])
    for i, H in enumerate(tgtnodes):
        add_gate(7.5 + i * 1e-6, FULL, [(tgt_of[(Q, H)], [sw[r] for r in byH[H]]) for Q in members[H] if tgt_of[(Q, H)][0] == 'b'])
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
