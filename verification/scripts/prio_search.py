import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import subprocess, random, sys, math
BIN = f"{ROOT}/csrc/tree"
h = int(sys.argv[1]); iters = int(sys.argv[2]); seed = int(sys.argv[3])
random.seed(seed)
def ev(pt):
    inp = f"{h}\n{' '.join(map(str, range(h)))}\n{' '.join(map(str, pt))}\n"
    return int(subprocess.run([BIN], input=inp, capture_output=True, text=True).stdout.split()[0])
cur = [h - 1 - x for x in range(h)]; cv = ev(cur); best = (cv, cur[:])
print("start", cv, flush=True)
T0 = cv * 0.002
for it in range(iters):
    nxt = cur[:]
    i, j = random.sample(range(h), 2)
    if random.random() < 0.5: nxt[i], nxt[j] = nxt[j], nxt[i]
    else:
        i, j = sorted((i, j)); x = nxt.pop(i); nxt.insert(j, x)
    nv = ev(nxt)
    T = T0 * (1 - it / iters) + 1e-9
    if nv <= cv or random.random() < math.exp((cv - nv) / T):
        cur, cv = nxt, nv
        if cv < best[0]: best = (cv, cur[:]); print(it, cv, cur, flush=True)
print("best", best)
