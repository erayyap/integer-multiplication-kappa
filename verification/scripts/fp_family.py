"""k-set designs: M_ST = A(|S cap T|) mod 2, labels Gram = B(|S cap T|) mod p (or over Q).
A, B arbitrary functions on 0..k given in binomial basis (any function on 0..k).
Need: A(k)=1, B(k)!=0, and for j<k: A(j)=1 => B(j)=0.
c = sum_{i: a_i=1} C(n,i);  xi <= sum_{i: b_i != 0} C(n,i)  (upper bound).
eta = (v - 6 c xi) / ((2v + 3 v z + 3c) xi^3);  theta = -ln(1-eta)/ln(xi^3)."""
from math import comb, log
import itertools, sys

def binom_mod(vals, p):
    b = []
    for j in range(len(vals)):
        s = sum(b[i] * comb(j, i) for i in range(j)) % p
        b.append((vals[j] - s) % p)
    return b

_cache = {}
def supports(Aset, k, p):
    key = (tuple(Aset), k, p)
    if key in _cache: return _cache[key]
    free = [j for j in range(k) if j not in Aset]
    sup = {}
    for fv in itertools.product(range(p), repeat=len(free)):
        vals = [0] * (k + 1)
        for j, x in zip(free, fv): vals[j] = x
        vals[k] = 1
        b = binom_mod(vals, p)
        sp = frozenset(i for i in range(k + 1) if b[i])
        sup.setdefault(sp, b)
    mins = [(sp, b) for sp, b in sup.items() if not any(o < sp for o in sup)]
    _cache[key] = mins
    return mins

def best_B(Aset, k, p, n):
    best = None
    for sp, b in supports(Aset, k, p):
        xi = sum(comb(n, i) for i in sp)
        if best is None or xi < best[0]:
            best = (xi, b)
    return best

def run(kmax, nmax, primes):
    res = []
    for k in range(2, kmax + 1):
        for bits in itertools.product((0, 1), repeat=k):
            Avals = list(bits) + [1]
            a = binom_mod(Avals, 2)
            Aset = [j for j in range(k) if Avals[j]]
            for p in primes:
                for n in range(2 * k, nmax + 1):
                    v = comb(n, k)
                    c = sum(comb(n, i) for i in range(k + 1) if a[i])
                    xb = best_B(Aset, k, p, n) if p * 0 == 0 else None
                    if xb is None: continue
                    xi = xb[0]
                    z = sum(comb(k, j) * comb(n - k, k - j) for j in Aset)
                    if v - 6 * c * xi <= 0 or xi < 2: continue
                    eta = (v - 6 * c * xi) / ((2 * v + 3 * v * z + 3 * c) * xi ** 3)
                    th = -log(1 - eta) / log(xi ** 3)
                    res.append((th, k, tuple(Avals), p, n, c, xi, v, z, tuple(xb[1])))
    res.sort(reverse=True)
    seen = set(); out = []
    for r in res:
        key = (r[1], r[2], r[3])
        if key in seen: continue
        seen.add(key); out.append(r)
    return out

if __name__ == "__main__":
    kmax = int(sys.argv[1]); nmax = int(sys.argv[2])
    primes = [int(x) for x in sys.argv[3].split(",")]
    for r in run(kmax, nmax, primes)[:20]:
        th, k, A, p, n, c, xi, v, z, b = r
        print(f"log2th={log(th,2):.3f} k={k} A={A} p={p} n={n} c={c} xi<={xi} v={v} z={z} B={b}")
