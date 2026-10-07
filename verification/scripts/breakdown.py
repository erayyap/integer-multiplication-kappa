import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import itertools, sys
from math import comb
sys.path.insert(0, f'{ROOT}/scripts')
from groupings import design
h = int(sys.argv[1]); k = int(sys.argv[2]); kind = sys.argv[3]
gS, gT = design(h, k, kind)
def info(groups):
    d = {}
    for gi, g in enumerate(groups):
        ker = set(g[0])
        for t in g[1:]: ker &= set(t)
        for t in g:
            pet = [x for x in t if x not in ker][0] if len(g) > 1 else None
            d[t] = (gi, len(g), pet)
    return d
iS, iT = info(gS), info(gT)
from collections import Counter
C = Counter(); sizes = Counter(len(g) for g in gS)
for X in itertools.combinations(range(h), 3):
    for S in itertools.combinations(range(h), 3):
        I = set(X) & set(S)
        if len(I) != 1: continue
        q = I.pop()
        a = 'P' if iS[X][2] == q else 'K'
        b = 'P' if iT[S][2] == q else 'K'
        C[a + b] += 1
tot = sum(C.values())
print({k_: round(v_/tot, 4) for k_, v_ in C.items()}, "group sizes", sorted(sizes.items()))
