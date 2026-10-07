"""F2 rank of the disjointness matrix on 2-subsets of [n] (lower bound on #rects per kernel,
odd covers allowed), and Q-rank (exact covers)."""
import itertools, sys
import numpy as np
def f2rank(M):
    M = M.copy() % 2; r = 0; rows, cols = M.shape
    for c in range(cols):
        piv = np.nonzero(M[r:, c])[0]
        if len(piv) == 0: continue
        pr = r + piv[0]; M[[r, pr]] = M[[pr, r]]
        nz = np.nonzero(M[:, c])[0]; nz = nz[nz != r]
        M[nz] ^= M[r]; r += 1
        if r == rows: break
    return r
for n in map(int, sys.argv[1:]):
    P = list(itertools.combinations(range(n), 2))
    D = np.array([[0 if set(a) & set(b) else 1 for b in P] for a in P], dtype=np.uint8)
    print(n, 'pairs', len(P), 'F2 rank', f2rank(D), 'Q rank', np.linalg.matrix_rank(D.astype(float)))
