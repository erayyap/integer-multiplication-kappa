"""theta_b and kappa ~ 0.0998 theta^2 if side wires per element were f * C(h,3) (per stage invocation)."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys
from math import comb, log
sys.path.insert(0, f'{ROOT}/scripts')
from networks import counts_bit_staged_side, theta
for f in (26.2, 13.1, 7.0, 3.0, 0.0):
    best = None
    for h1 in range(12, 80):
        for h2 in range(12, 80):
            st = counts_bit_staged_side((h1, h2, h1), (round(f * comb(h1, 3)), round(f * comb(h2, 3)), round(f * comb(h1, 3))))
            if 2 * st['L'] >= st['N']: continue
            t = theta(st['eta_b'], st['m'])
            if best is None or t > best[0]: best = (t, h1, h2)
    t, h1, h2 = best
    print(f"side/v={f:5.1f}  best h=({h1},{h2},{h1})  log2 theta={log(t,2):.3f}  log2 kappa~{log(0.0998*t*t,2):.3f}")
