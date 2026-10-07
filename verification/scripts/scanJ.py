import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys, ast
from math import comb, log
sys.path.insert(0, f'{ROOT}/scripts')
from networks import counts_bit_staged_side, theta
def cross(h):  # rot order: total src+tgt chains = 2h(h-3)^2 (exact h=9..19), per side h(h-3)^2
    return h * (h - 3) ** 2
def cross_side(h):
    d = ast.literal_eval(open(f'{ROOT}/data/sorted_{h-1}.txt').read().strip().splitlines()[-1])
    return h * d['rects'] + 2 * cross(h) - 2 * comb(h, 3)
res = []
for h1 in range(31, 58):
    for h2 in range(31, 58):
        st = counts_bit_staged_side((h1, h2, h1), (cross_side(h1), cross_side(h2), cross_side(h1)))
        if 2 * st['L'] >= st['N']: continue
        res.append((theta(st['eta_b'], st['m']), h1, h2))
res.sort(reverse=True)
for t, h1, h2 in res[:5]: print(h1, h2, 'log2 theta', round(log(t, 2), 4))
