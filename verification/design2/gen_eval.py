"""Generic evaluator for a 3-stage bit network built from a per-coordinate design.

Per coordinate j: v_j elements, label dim f_j (1-dim labels t_S, k=1), centre F2-rank c_j,
side_j = side wires (+copy/accumulator chains) per invocation.  Stage 1->3 scratch reuse (needs
equal stage-1/3 designs).  L = sum inv_j c_j f_j (centre wire loses P(x)F(x)Q, dim f_j).
theta = -ln(1-eta)/ln m.  Returns log2 theta.
"""
from fractions import Fraction as Fr
from math import log, comb

def theta(designs, reuse=True):
    vs = [d['v'] for d in designs]; fs = [d['f'] for d in designs]
    N = vs[0] * vs[1] * vs[2]; m = fs[0] * fs[1] * fs[2]
    inv = [N // vs[j] for j in range(3)]
    scratch = [inv[j] * (designs[j]['side'] + designs[j]['c']) for j in range(3)]
    W = 2 * N + (scratch[0] + scratch[1] if reuse else sum(scratch))
    L = sum(inv[j] * designs[j]['c'] * fs[j] for j in range(3))
    s = W * m - N + 2 * L
    eta = (W * m - s) / (W * m)
    if eta <= 0: return None
    return log(-__import__("math").log1p(-eta) / log(m), 2)

def sets_design(h, s, cdeg, side_per_v):
    v = comb(h, s)
    return dict(v=v, f=h, c=comb(h, cdeg), side=side_per_v * v)

if __name__ == "__main__":
    import sys
    # reference: triples at K: side per invocation from verify (side/v = 20.565 at h=51), c = h
    def tri(h, sv): v = comb(h, 3); return dict(v=v, f=h, c=h, side=sv * v)
    print('triples K-like (51,50,51) side/v 20.565:', theta([tri(51, 20.565), tri(50, 20.565), tri(51, 20.565)]))
    target = theta([tri(51, 20.565), tri(50, 20.565), tri(51, 20.565)])
    # 7-sets, centre C(|S cap T|,3) mod 2 (rank <= C(h,3)), labels f = h
    for s, cd in ((7, 3), (15, 7)):
        print(f'--- {s}-sets, centre degree {cd}')
        for h in range(14, 61, 2):
            # break-even side/v
            lo, hi = 0.0, 1e7
            if theta([sets_design(h, s, cd, 0)] * 3) is None or theta([sets_design(h, s, cd, 0)] * 3) < target:
                print(h, 'never beats', theta([sets_design(h, s, cd, 0)] * 3)); continue
            for _ in range(60):
                mid = (lo + hi) / 2
                t = theta([sets_design(h, s, cd, mid)] * 3)
                if t is not None and t > target: lo = mid
                else: hi = mid
            print(h, 'theta(side=0)=%.3f' % theta([sets_design(h, s, cd, 0)] * 3), 'theta(side=20v)=%.3f' % theta([sets_design(h, s, cd, 20)] * 3), 'break-even side/v = %.1f' % lo)
