"""Search small 'stage designs' for the bit motif network.

A design: elements S with (x_S, y_S) in F_2^c (x_S . y_S = 1) and a label t_S in
Q^h (non-isotropic for a fixed nondegenerate integer form Q).  M_ST = x_S . y_T
(mod 2).  Required: S != T and M_ST = 1  =>  <t_S, t_T>_Q = 0.
Compatible pairs form a graph; a design is a clique.  We report max v and
the ratio R = v / (c h).
"""
import itertools, sys
import numpy as np
import networkx as nx


def primitive(v):
    from math import gcd
    g = 0
    for a in v:
        g = gcd(g, abs(a))
    if g == 0:
        return None
    v = tuple(a // g for a in v)
    for a in v:
        if a != 0:
            if a < 0:
                v = tuple(-b for b in v)
            break
    return v


def labels(h, B, Q):
    out = set()
    for v in itertools.product(range(-B, B + 1), repeat=h):
        p = primitive(v)
        if p is None:
            continue
        a = np.array(p)
        if a @ Q @ a != 0:
            out.add(p)
    return sorted(out)


def types(c):
    vecs = [v for v in itertools.product((0, 1), repeat=c) if any(v)]
    return [(x, y) for x in vecs for y in vecs if sum(a * b for a, b in zip(x, y)) % 2 == 1]


def dot2(x, y):
    return sum(a * b for a, b in zip(x, y)) % 2


def build(c, h, B, Q):
    L = labels(h, B, Q)
    T = types(c)
    La = np.array(L)
    G = La @ Q @ La.T
    nodes = [(ti, li) for ti in range(len(T)) for li in range(len(L))]
    g = nx.Graph()
    g.add_nodes_from(range(len(nodes)))
    for i in range(len(nodes)):
        ti, li = nodes[i]
        for j in range(i + 1, len(nodes)):
            tj, lj = nodes[j]
            orth = G[li, lj] == 0
            if orth or (dot2(T[ti][0], T[tj][1]) == 0 and dot2(T[tj][0], T[ti][1]) == 0):
                g.add_edge(i, j)
    return g, nodes, T, L


if __name__ == "__main__":
    c, h, B = map(int, sys.argv[1:4])
    form = sys.argv[4] if len(sys.argv) > 4 else "I"
    if form == "I":
        Q = np.eye(h, dtype=int)
    elif form == "L":
        Q = np.diag([1] * (h - 1) + [-1])
    g, nodes, T, L = build(c, h, B, Q)
    print("nodes", g.number_of_nodes(), "edges", g.number_of_edges())
    cl, w = nx.max_weight_clique(g, weight=None)
    print("c,h,form", c, h, form, "max v", w, "R=v/(ch)", w / (c * h))
    E = 0
    for a in cl:
        for b in cl:
            if a != b and dot2(T[nodes[a][0]][0], T[nodes[b][0]][1]) == 1:
                E += 1
    print("|E|", E, "z", E / w)
    for a in cl:
        print(T[nodes[a][0]], L[nodes[a][1]])
