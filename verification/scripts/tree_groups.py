"""Laminar (binary-tree) side-wire groups.  Python reference + C counter (csrc/tree)."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import itertools, subprocess, sys
from math import comb
BIN = f"{ROOT}/csrc/tree"


def families(h, prio):
    fam = {}
    for a, b in itertools.combinations(range(h), 2):
        P = sorted(p for p in range(h) if p not in (a, b) and prio[p] < prio[a] and prio[p] < prio[b])
        if P: fam[(a, b)] = P
    return fam


def tree(n, lo=0, hi=None):
    hi = n if hi is None else hi
    node = (lo, hi)
    if hi - lo == 1: return (node, None, None)
    mid = (lo + hi) // 2
    return (node, tree(n, lo, mid), tree(n, mid, hi))


def cover(K, P, L, Q):
    """returns list of rectangles ((lo,hi),(lo,hi)) in petal-index space."""
    nb = [[len((set(K) | {p}) & (set(L) | {q})) == 1 for q in Q] for p in P]
    memo = {}
    def f(A, B):
        key = (A[0], B[0])
        if key in memo: return memo[key]
        (a0, a1), (b0, b1) = A[0], B[0]
        c = sum(nb[x][y] for x in range(a0, a1) for y in range(b0, b1))
        if c == 0: r = (0, [])
        elif c == (a1 - a0) * (b1 - b0): r = (1, [(A[0], B[0])])
        else:
            opts = []
            if A[1]:
                x, y = f(A[1], B), f(A[2], B); opts.append((x[0] + y[0], x[1] + y[1]))
            if B[1]:
                x, y = f(A, B[1]), f(A, B[2]); opts.append((x[0] + y[0], x[1] + y[1]))
            r = min(opts, key=lambda t: t[0])
        memo[key] = r
        return r
    return f(tree(len(P)), tree(len(Q)))


def count(h, pS, pT):
    inp = f"{h}\n{' '.join(map(str, pS))}\n{' '.join(map(str, pT))}\n"
    out = subprocess.run([BIN], input=inp, capture_output=True, text=True).stdout.split()
    return int(out[0]), int(out[1])


def anti(h):
    return list(range(h)), [h - 1 - x for x in range(h)]


if __name__ == "__main__":
    for h in map(int, sys.argv[1].split(',')):
        tot, pairs = count(h, *anti(h))
        v = comb(h, 3); z = 3 * comb(h - 3, 2)
        assert pairs == v * z
        print(h, tot, round(tot / (v * z), 4), flush=True)
