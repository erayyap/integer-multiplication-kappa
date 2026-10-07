"""Frankl-Wilson style stage designs for the bit motif network.

Elements: k-subsets of [n].  F_2 coefficient M_ST = f(|S cap T|) (f(k)=1),
realized with c centers, c = sum_{i: a_i=1} C(n,i) where f(j)=sum a_i C(j,i) mod 2.
Labels: Gram g(|S cap T|) over Q, g = prod_{a in Z}(j-a), Z = side intersection sizes.
h = rank_Q of the Gram matrix via Johnson-scheme eigenvalues.
eta = (v - 6 c h) / ((2v + 3 v z + 3c) h^3);  theta = -ln(1-eta)/ln(h^3).
"""
from math import comb, log
from fractions import Fraction as Fr
import itertools, sys


def binom_coeffs_mod2(fvals):
    # f(j) = sum_i a_i C(j,i) mod 2, solve triangular system
    k = len(fvals) - 1
    a = []
    for j in range(k + 1):
        s = sum(a[i] * comb(j, i) for i in range(j)) % 2
        a.append((fvals[j] - s) % 2)
    return a


def poly_to_binom(coeffs_j):
    """polynomial in j (list of power coeffs) -> coefficients b_i with g(j)=sum b_i C(j,i)."""
    deg = len(coeffs_j) - 1
    vals = [sum(Fr(cf) * j ** p for p, cf in enumerate(coeffs_j)) for j in range(deg + 1)]
    b = []
    for j in range(deg + 1):
        s = sum(b[i] * comb(j, i) for i in range(j))
        b.append(vals[j] - s)
    return b


def rank_gram(n, k, b):
    """rank over Q of sum_i b_i W_i^T W_i on k-subsets of [n] (k <= n/2)."""
    h = 0
    for r in range(k + 1):
        lam = sum(b[i] * comb(k - r, i - r) * comb(n - i - r, k - i)
                  for i in range(r, len(b)))
        if lam != 0:
            h += comb(n, r) - (comb(n, r - 1) if r > 0 else 0)
    return h


def evaluate(n, k, fvals):
    assert fvals[k] == 1
    Z = [j for j in range(k) if fvals[j] == 1]
    a = binom_coeffs_mod2(fvals)
    c = sum(comb(n, i) for i in range(k + 1) if a[i])
    # g(j) = prod (j - z)
    poly = [Fr(1)]
    for z0 in Z:
        new = [Fr(0)] * (len(poly) + 1)
        for p, cf in enumerate(poly):
            new[p + 1] += cf
            new[p] -= z0 * cf
        poly = new
    b = poly_to_binom(poly)
    h = rank_gram(n, k, b)
    v = comb(n, k)
    zz = sum(comb(k, j) * comb(n - k, k - j) for j in Z)
    if v - 6 * c * h <= 0:
        return None
    eta = Fr(v - 6 * c * h, (2 * v + 3 * v * zz + 3 * c) * h ** 3)
    m = h ** 3
    if m < 2:
        return None
    th = -log(1 - float(eta)) / log(m)
    return dict(n=n, k=k, f=fvals, Z=Z, c=c, h=h, v=v, z=zz, eta=eta, theta=th)


if __name__ == "__main__":
    # sanity: paper triples (form I - J/9 has rank n)
    r = evaluate(100, 3, [0, 1, 0, 1])
    print("triples n=100:", r["c"], r["h"], float(r["eta"]), log(r["theta"], 2))
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 120
    best = []
    for k in range(1, kmax + 1):
        for bits in itertools.product((0, 1), repeat=k):
            fvals = list(bits) + [1]
            for n in range(2 * k, nmax + 1):
                r = evaluate(n, k, fvals)
                if r:
                    best.append((r["theta"], r))
    best.sort(key=lambda t: -t[0])
    seen = set()
    for th, r in best:
        key = (r["k"], tuple(r["f"]))
        if key in seen:
            continue
        seen.add(key)
        print(f"log2theta={log(th,2):.3f} n={r['n']} k={r['k']} f={r['f']} Z={r['Z']} c={r['c']} h={r['h']} v={r['v']} z={r['z']}")
        if len(seen) > 15:
            break
