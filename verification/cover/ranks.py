import sys, itertools, numpy as np
import galois
n = int(sys.argv[1]); GF = galois.GF(2)
pairs = list(itertools.combinations(range(n), 2)); K = len(pairs)
def mat(f):
    M = np.zeros((K, K), dtype=np.uint8)
    for i, (a, b) in enumerate(pairs):
        for j, (c, d) in enumerate(pairs):
            if len({a, b, c, d}) == 4 and f(a, b, c, d): M[i, j] = 1
    return GF(M)
R = {'D': lambda a,b,c,d: True,
     'R1 p1<q1,p2<q2': lambda a,b,c,d: a<c and b<d,
     'cross PQPQ': lambda a,b,c,d: a<c<b<d,
     'sep PPQQ': lambda a,b,c,d: b<c,
     'R3 PQQP': lambda a,b,c,d: a<c<d<b,
     'R1+R3 (p1<q1)': lambda a,b,c,d: a<c,
     'cross+nest (interleaved)': lambda a,b,c,d: not (b<c or d<a)}
for k, f in R.items(): print(n, k, np.linalg.matrix_rank(mat(f)), 'C(n,2)=', K, flush=True)
