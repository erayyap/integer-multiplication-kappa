"""Exact counts for the paper's two motif networks (Section 3) as functions of h."""
from fractions import Fraction as Fr
from math import comb, log


def counts(h):
    v = comb(h, 3)
    N = v ** 3
    m = h ** 3
    I = 3 * v ** 2
    zb = 3 * comb(h - 3, 2)
    zc = comb(h - 3, 3) + 3 * (h - 3)
    cb, cc = h, h + 1
    Wb = 2 * N + I * (v * zb + cb)
    Wc = 2 * N + I * (v * zc + cc)
    Lb = I * cb * h
    Lc = I * cc * h
    sb = Wb * m - N + 2 * Lb
    sc = Wc * m - 2 * N + 2 * Lc
    return dict(h=h, v=v, N=N, m=m, I=I, Wb=Wb, Wc=Wc, sb=sb, sc=sc, Lb=Lb, Lc=Lc,
                eta_b=Fr(Wb * m - sb, Wb * m), eta_c=Fr(Wc * m - sc, Wc * m))


def theta(eta, m):
    """sup of 1-tau with s/W < m^tau, i.e. -ln(1-eta)/ln m (float)."""
    return -log(1 - float(eta)) / log(m)


if __name__ == "__main__":
    c = counts(100)
    assert c["Wb"] == 177176569091445000000 and c["Wc"] == 1873807244643542670000
    assert c["sb"] == 177176569088785861287000000 and c["sc"] == 1873807244636671267308000000
    assert c["eta_b"] == Fr(339, 22587335000000) and c["eta_c"] == Fr(73, 19906842167500)
    print("paper h=100 reproduced")
    best = {}
    for h in range(7, 200):
        c = counts(h)
        if h == 9:
            continue
        for k in ("b", "c"):
            e = c["eta_" + k]
            if e <= 0:
                continue
            t = theta(e, c["m"])
            if t > best.get(k, (0,))[0]:
                best[k] = (t, h)
    for k, (t, h) in best.items():
        print(k, "best h", h, "theta=1-tau", t, "log2", log(t, 2))
    c = counts(100)
    print("h=100: theta_b", log(theta(c["eta_b"], c["m"]), 2), "theta_c", log(theta(c["eta_c"], c["m"]), 2))


def counts_bit_staged(hs, reuse=True):
    """Bit network with coordinate j using the triple design on hs[j] points.

    reuse=True: every stage-3 scratch wire is a stage-1 scratch wire continued
    (requires hs[0]==hs[2]); its stage-1 final label F(x)t_{a2}(x)t_{a3} lies in
    its stage-3 initial label (t_{a1'}(x)t_{a2'})^perp (x) <t_{a3}> whenever
    a2, a2' are neighbours, so the continuation edge is increasing (no loss).
    """
    h1, h2, h3 = hs
    vs = [comb(h, 3) for h in hs]
    zs = [3 * comb(h - 3, 2) for h in hs]
    cs = list(hs)
    N = vs[0] * vs[1] * vs[2]
    m = h1 * h2 * h3
    inv = [N // vs[j] for j in range(3)]
    scratch = [inv[j] * (vs[j] * zs[j] + cs[j]) for j in range(3)]
    if reuse:
        assert h1 == h3 and scratch[0] == scratch[2]
        W = 2 * N + scratch[0] + scratch[1]
    else:
        W = 2 * N + sum(scratch)
    L = sum(inv[j] * cs[j] * hs[j] for j in range(3))
    s = W * m - N + 2 * L
    # caller checks 2L < N (needed for a saving)
    return dict(N=N, m=m, W=W, L=L, s=s, eta_b=Fr(W * m - s, W * m))


def counts_bit_staged_side(hs, sides, reuse=True):
    """As counts_bit_staged, but stage j uses sides[j] side wires per invocation
    (grouped side-wire gates, scripts/groupings.py) instead of v_j*z_j."""
    h1, h2, h3 = hs
    vs = [comb(h, 3) for h in hs]
    cs = list(hs)
    N = vs[0] * vs[1] * vs[2]
    m = h1 * h2 * h3
    inv = [N // vs[j] for j in range(3)]
    scratch = [inv[j] * (sides[j] + cs[j]) for j in range(3)]
    if reuse:
        assert h1 == h3 and scratch[0] == scratch[2]
        W = 2 * N + scratch[0] + scratch[1]
    else:
        W = 2 * N + sum(scratch)
    L = sum(inv[j] * cs[j] * hs[j] for j in range(3))
    s = W * m - N + 2 * L
    return dict(N=N, m=m, W=W, L=L, s=s, eta_b=Fr(W * m - s, W * m))
