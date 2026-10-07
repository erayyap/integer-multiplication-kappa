import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys, math, glob
sys.path.insert(0, f'{ROOT}/scripts')
from networks import counts_bit_staged_side, theta
side = {}
for f in glob.glob(f'{ROOT}/data/tree_*.txt'):
    h = int(f.split('_')[-1].split('.')[0]); side[h] = int(open(f).read().split()[0])
best = None; rows = []
for h1 in sorted(side):
    for h2 in sorted(side):
        st = counts_bit_staged_side((h1, h2, h1), (side[h1], side[h2], side[h1]))
        if 2 * st['L'] >= st['N']: continue
        t = theta(st['eta_b'], st['m'])
        rows.append((t, h1, h2))
rows.sort(reverse=True)
for t, h1, h2 in rows[:8]: print(round(math.log2(t), 4), h1, h2)
