import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys, random, math, time
sys.path.insert(0, f'{ROOT}/scripts')
from cross_var import run
h = int(sys.argv[1]); seed = int(sys.argv[2]); T0 = float(sys.argv[3]); random.seed(seed)
ords = {q: [(q + 1 + i) % h for i in range(h - 1)] for q in range(h)}
def ev(o): return run(h, lambda q, hh: o[q])
cur = ev(ords); best = cur; print('start', cur, flush=True)
t0 = time.time(); it = 0
while time.time() - t0 < 1500:
    it += 1; T = T0 * (1 - (time.time() - t0) / 1500) + 1e-3
    q = random.randrange(h); o = ords[q][:]; i, j = random.sample(range(h - 1), 2)
    if random.random() < 0.5: o[i], o[j] = o[j], o[i]
    else: x = o.pop(i); o.insert(j, x)
    new = dict(ords); new[q] = o; v = ev(new)
    if v <= cur or random.random() < math.exp((cur - v) / T):
        ords, cur = new, v
        if cur < best: best = cur; print(it, 'best', best, {k: ords[k] for k in range(h)}, flush=True)
print('done its', it, 'best', best, flush=True)
