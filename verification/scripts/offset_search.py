import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys, random, time
sys.path.insert(0, f'{ROOT}/scripts')
from cross_var import run
h = int(sys.argv[1]); seed = int(sys.argv[2]); random.seed(seed)
def ev(d): return run(h, lambda q, hh: [(q + x) % hh for x in d])
d = list(range(1, h)); best = ev(d); print('start rot', best, flush=True)
t0 = time.time(); it = 0
while time.time() - t0 < 4800:
    e = d[:]; i, j = random.sample(range(h - 1), 2)
    if random.random() < 0.5: e[i], e[j] = e[j], e[i]
    else: x = e.pop(i); e.insert(j, x)
    v = ev(e); it += 1
    if v <= best:
        if v < best: print(it, 'improved', v, e, flush=True)
        best, d = v, e
print('final', best, d, flush=True)
