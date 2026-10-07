import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys, math
from math import comb
sys.path.insert(0, f'{ROOT}/scripts')
from groupings import design, count
from networks import counts_bit_staged_side, counts_bit_staged, theta
side = {}
def best_side(h, ks=(3,4,5,6,7)):
    if h not in side:
        side[h] = min((count(h, *design(h, k, "anti"))[0], k) for k in ks)
    return side[h]
best = None
for h1 in range(36, 62, 2):
    for h2 in range(30, 62, 2):
        s1, s2 = best_side(h1), best_side(h2)
        st = counts_bit_staged_side((h1, h2, h1), (s1[0], s2[0], s1[0]))
        if 2 * st['L'] >= st['N']: continue
        t = theta(st['eta_b'], st['m'])
        if best is None or t > best[0]: best = (t, h1, h2, s1, s2)
    print(h1, best and (math.log2(best[0]), best[1:]), flush=True)
st = counts_bit_staged((48, 42, 48))
print("E baseline", math.log2(theta(st['eta_b'], st['m'])))
