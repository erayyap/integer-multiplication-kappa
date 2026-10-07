"""Side-wire compression by sunflower groups (bit network, one stage, triples on [h]).

Each triple S gets a kernel pair; triples with the same kernel are chunked into
groups of <= k petals.  Gates: source X has own gate P(x)t_X and group gate
P(x)span(t_G'); target S has group gate P(x)(span t_G)^perp and own gate
P(x)t_S^perp.  A wire may be (own|group source) x (own|group target) provided the
source label is contained in the target label, i.e. every source in it is a
neighbour (|X cap S| = 1) of every target in it.  Each neighbour pair must be
covered exactly once.  Per block (G',G) we pick the cheapest of:
complete -> 1 wire; rows -> full-row wires + singletons; cols -> symmetric;
all singletons.  Returns number of side wires per invocation."""
import itertools, sys
import numpy as np
from math import comb


def triples(h):
    T = list(itertools.combinations(range(h), 3))
    return T


def assign(h, k, mode="sum3"):
    T = triples(h)
    bykernel = {}
    for t in T:
        a, b, c = t
        if mode == "sum3":
            drop = t[(a + b + c) % 3]
        elif mode == "max":
            drop = c
        ker = tuple(x for x in t if x != drop)
        bykernel.setdefault(ker, []).append(drop)
    groups = []
    for ker, pet in bykernel.items():
        pet.sort()
        for i in range(0, len(pet), k):
            groups.append([tuple(sorted(ker + (p,))) for p in pet[i:i + k]])
    return T, groups


def count(h, k, mode="sum3", verbose=False):
    T, groups = assign(h, k, mode)
    v = len(T)
    idx = {t: i for i, t in enumerate(T)}
    mask = np.array([(1 << a) | (1 << b) | (1 << c) for a, b, c in T], dtype=np.int64)
    # order sources by group
    order = [idx[t] for g in groups for t in g]
    gid = np.concatenate([[gi] * len(g) for gi, g in enumerate(groups)])
    gsize = np.array([len(g) for g in groups])
    G = len(groups)
    smask = mask[order]
    pc = np.vectorize(lambda x: bin(x).count("1"))
    # popcount table for 3-bit intersections: use bit tricks via numpy
    total = 0
    singles = 0
    starts = np.concatenate([[0], np.cumsum(gsize)])
    for gt in range(G):
        tg = [mask[idx[t]] for t in groups[gt]]
        kk = len(tg)
        nb = np.zeros((v, kk), dtype=bool)
        for j, m in enumerate(tg):
            inter = smask & m
            # popcount of inter (at most 3 bits set): count via x&(x-1)
            x = inter
            c = (x != 0).astype(np.int8)
            x = x & (x - 1)
            c += (x != 0)
            x = x & (x - 1)
            c += (x != 0)
            nb[:, j] = (c == 1)
        # per source group block
        rows_full = nb.all(axis=1)
        npairs_row = nb.sum(axis=1)
        for gs in range(G):
            a0, a1 = starts[gs], starts[gs + 1]
            blk = nb[a0:a1]
            n = int(npairs_row[a0:a1].sum())
            if n == 0:
                continue
            if n == blk.size:
                total += 1
                continue
            fr = rows_full[a0:a1]
            cost_rows = int(fr.sum()) + int(npairs_row[a0:a1][~fr].sum())
            fc = blk.all(axis=0)
            cost_cols = int(fc.sum()) + int(blk[:, ~fc].sum())
            total += min(n, cost_rows, cost_cols)
    z = 3 * comb(h - 3, 2)
    return v, total, v * z


if __name__ == "__main__":
    h = int(sys.argv[1]); k = int(sys.argv[2]); mode = sys.argv[3] if len(sys.argv) > 3 else "sum3"
    v, tot, base = count(h, k, mode)
    print(f"h={h} k={k} mode={mode} v={v} wires={tot} baseline={base} ratio={tot/base:.4f} z_eff={tot/v:.1f}")
