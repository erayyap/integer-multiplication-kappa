import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import ast, sys
from math import comb, log
sys.path.insert(0, f'{ROOT}/scripts')
from networks import counts_bit_staged_side, theta
def copy_side(h):
    d = ast.literal_eval(open(f'{ROOT}/data/copy_{h-1}.txt').read().strip().splitlines()[-1])
    return h * d['total'] - 2 * comb(h, 3)
if __name__ == "__main__":
    res = []
    for h1 in range(37, 58):
        for h2 in range(37, 58):
            st = counts_bit_staged_side((h1, h2, h1), (copy_side(h1), copy_side(h2), copy_side(h1)))
            if 2 * st['L'] >= st['N']: continue
            res.append((log(theta(st['eta_b'], st['m']), 2), h1, h2))
    res.sort(reverse=True); print(res[:8])
    for h in (43, 48, 49): print(h, copy_side(h), copy_side(h) / comb(h, 3))
