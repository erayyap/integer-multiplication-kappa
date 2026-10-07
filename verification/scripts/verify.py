"""Exact rational verification of the improved parameter packages.

Every inequality below is checked with Fractions.  Logarithms enter only
through rigorous interval bounds (mpmath.iv) on ln m and log_m B.

Packages
  A3 : paper's algorithm unchanged; bit network h=46, complex network h=25,
       all fixed rationals re-chosen.  Cost table exactly as in the paper.
  B  : A3 + CRT reversal and Gaussian axis exposure by O(d) non-adjacent
       chunk interchanges (row 'CRT and axis layouts' becomes d(1+l^tau)).
       Guard exponent C1 = nu+3 = 14, as stated in the paper (eps*C1<1).
  C  : B + delta small (any fixed 0<delta<1/8 is allowed by Lemma 7.x),
       the convenience bound eps<1/12 relaxed to eps<1/4 (gamma=o(b)),
       and the guard depth bound taken with the leaf size d^beta:
       Delta = d^(beta+(1-beta)log_m B + o(1)), requiring
       eps*(beta+(1-beta)*log_m B) < 1.
"""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
from fractions import Fraction as F
from math import comb
import mpmath
from mpmath import iv
iv.dps = 60
import sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from networks import counts, counts_bit_staged, counts_bit_staged_side
from groupings import design, count as side_count
import tree_groups


def _ep(t):
    sign, man, exp, bc = t
    v = F(man) * (F(2)**exp if exp >= 0 else F(1, 2**(-exp)))
    return -v if sign else v


def up(I):
    return _ep(I._mpi_[1])


def lo(I):
    return _ep(I._mpi_[0])


def ln_upper(x):
    return up(iv.log(iv.mpf(x)))


def exponent_ok(theta, eta, m):
    """sufficient: theta*ln m < eta  ==>  m^(1-theta) > m(1-theta ln m) > m(1-eta) = s/W."""
    U = ln_upper(m)
    # also need m^{-theta} > 1 - theta ln m, true since e^{-x} > 1-x.
    return theta * U < eta


def dyadic_below(x, bits=12):
    """largest k*2^-(e+bits) <= x with 2^-e ~ x, as Fraction."""
    x = F(x)
    e = 0
    while F(1, 2**e) > x:
        e += 1
    den = 2**(e + bits)
    return F(int(x * den) - 1, den)


BIT_STAGES = (48, 42, 48)
# package F: grouped side-wire gates (anti sunflower groupings, group size k)
F_STAGES = ((48, 5), (44, 5), (48, 5))


# package G: laminar binary-tree side-wire groups per kernel (csrc/tree, scripts/tree_groups.py)
G_STAGES = (49, 43, 49)


def tree_side_wires(h):
    tot, pairs = tree_groups.count(h, *tree_groups.anti(h))
    assert pairs == comb(h, 3) * 3 * comb(h - 3, 2)
    return tot


# package H: copy-wire side covers per kernel point (scripts/copy_cover.py, check_copy.py)
H_STAGES = (51, 48, 51)


def copy_side_wires(h):
    import ast
    d = ast.literal_eval(open(f'{ROOT}/data/copy_{h-1}.txt').read().strip().splitlines()[-1])
    assert d['n'] == h - 1 and d['cells'] == comb(h - 1, 2) * comb(h - 3, 2)
    return h * d['total'] - 2 * comb(h, 3)


# package I: copy wires with the sorted-order cover (scripts/copy_sorted.py, check_copy.py)
I_STAGES = (51, 50, 51)


def sorted_side_wires(h):
    import ast
    d = ast.literal_eval(open(f'{ROOT}/data/sorted_{h-1}.txt').read().strip().splitlines()[-1])
    assert d['n'] == h - 1 and d['cells'] == comb(h - 1, 2) * comb(h - 3, 2)
    return h * d['total'] - 2 * comb(h, 3)


# package J: package I cover, copy/accumulator chains across the 3 kernels of a triple
# (scripts/cross_width.py, check_cross.py).  cross(h) = per-invocation source chains,
# exact for h = 9..34 (data/cross_width*.log); targets identical.
J_STAGES = (51, 50, 51)


def cross_chains(h):
    v = 7 * h ** 3 - 39 * h ** 2 + 26 * h + 102
    assert v % 6 == 0
    return v // 6


def cross_side_wires(h):
    import ast
    d = ast.literal_eval(open(f'{ROOT}/data/sorted_{h-1}.txt').read().strip().splitlines()[-1])
    assert d['n'] == h - 1 and d['cells'] == comb(h - 1, 2) * comb(h - 3, 2)
    return h * d['rects'] + 2 * cross_chains(h) - 2 * comb(h, 3)


# package K: as J but kernel q orders the other points cyclically q+1, q+2, ... (mod h);
# total cross-kernel chains (source + target) = 2 h (h-3)^2, exact for h = 9..19 (data/cross_rot_*.log);
# explicit networks: check_cross.py with ORDER=rot.
K_STAGES = (51, 50, 51)


def rot_side_wires(h):
    import ast
    d = ast.literal_eval(open(f'{ROOT}/data/sorted_{h-1}.txt').read().strip().splitlines()[-1])
    assert d['n'] == h - 1 and d['cells'] == comb(h - 1, 2) * comb(h - 3, 2)
    return h * d['rects'] + 2 * h * (h - 3) ** 2 - 2 * comb(h, 3)


def side_wires(h, k):
    tot, pairs, _ = side_count(h, *design(h, k, "anti"))
    assert pairs == comb(h, 3) * 3 * comb(h - 3, 2)
    return tot


def check(pkg, bit_h=46, cx_h=25):
    nb, nc = counts(bit_h), counts(cx_h)
    eb, ec, mb, mc = nb['eta_b'], nc['eta_c'], nb['m'], nc['m']
    if pkg in ('D', 'E'):
        st = counts_bit_staged(BIT_STAGES, reuse=True)
        assert 2 * st['L'] < st['N'] and st['s'] < st['W'] * st['m']
        eb, mb = st['eta_b'], st['m']
    if pkg in ('F', 'G', 'H', 'I', 'J', 'K'):
        if pkg == 'K':
            st = counts_bit_staged_side(K_STAGES, tuple(rot_side_wires(h) for h in K_STAGES), reuse=True)
        elif pkg == 'J':
            st = counts_bit_staged_side(J_STAGES, tuple(cross_side_wires(h) for h in J_STAGES), reuse=True)
        elif pkg == 'I':
            st = counts_bit_staged_side(I_STAGES, tuple(sorted_side_wires(h) for h in I_STAGES), reuse=True)
        elif pkg == 'H':
            st = counts_bit_staged_side(H_STAGES, tuple(copy_side_wires(h) for h in H_STAGES), reuse=True)
        elif pkg == 'F':
            st = counts_bit_staged_side(tuple(h for h, _ in F_STAGES),
                                        tuple(side_wires(h, k) for h, k in F_STAGES), reuse=True)
        else:
            st = counts_bit_staged_side(G_STAGES, tuple(tree_side_wires(h) for h in G_STAGES), reuse=True)
        assert 2 * st['L'] < st['N'] and st['s'] < st['W'] * st['m']
        eb, mb = st['eta_b'], st['m']
    th_b = dyadic_below(eb / ln_upper(mb), 20)
    th_c = dyadic_below(ec / ln_upper(mc), 20)
    assert exponent_ok(th_b, eb, mb) and exponent_ok(th_c, ec, mc)
    tau, sigma = 1 - th_b, 1 - th_c
    # complex-network guard quantities
    Bc = nc['sc'] + 64 * (nc['Wc'] + mc + 1)**3
    nu = 1
    while mc**nu < Bc:
        nu += 1
    logmB_up = up(iv.log(iv.mpf(Bc)) / iv.log(iv.mpf(mc)))

    if pkg == 'A3':
        # g4 = th - eps(2-tau) forces eps ~ th; balance g2=eps c th, g3=eps(th - c/beta), g4
        beta = 1 - dyadic_below(F(1, 10**7), 8)
        c = dyadic_below(th_b * beta / (1 + th_b * beta) * (1 - th_b), 20)
        eps = dyadic_below(th_b / 2, 20)
        delta, C1 = F(1, 16), 20
    elif pkg == 'B':
        beta = 1 - dyadic_below(F(1, 10**7), 8)
        c = dyadic_below(th_b * beta / (1 + th_b * beta) * (1 - th_b), 20)
        C1 = nu + 3
        eps = F(1, C1) - F(1, 2**20)
        delta = F(1, 16)
    elif pkg in ('E', 'F', 'G', 'H', 'I', 'J', 'K'):
        beta = 1 - dyadic_below(F(1, 10**7), 8)
        c = dyadic_below(th_b * beta / (1 + th_b * beta) * (1 - th_b), 20)
        delta = F(1, 2**20)
        eps = F(1, 5) - F(1, 2**17)
        C1 = None
    elif pkg in ('C', 'D'):
        beta = 1 - dyadic_below(F(1, 10**7), 8)
        c = dyadic_below(th_b * beta / (1 + th_b * beta) * (1 - th_b), 20)
        delta = F(1, 2**20)
        eps = F(1, 6) - F(1, 2**17)
        C1 = None
    lam_lo = max(tau, sigma, tau * (1 + c / beta))
    lamp_lo = max(sigma + beta * (1 - sigma), lam_lo)
    gap = 1 - lamp_lo
    assert gap > 0
    lam = lam_lo + gap / 2**30
    lamp = max(lam, sigma + beta * (1 - sigma)) + gap / 2**30
    cons = {
        '0<tau<1': 0 < tau < 1,
        '0<sigma<1': 0 < sigma < 1,
        '0<beta<1': 0 < beta < 1,
        'max(tau,sigma)<lam<1': max(tau, sigma) < lam < 1,
        'tau(1+c/beta)<lam': tau * (1 + c / beta) < lam,
        'max(lam,sigma+beta(1-sigma))<lamp<1': max(lam, sigma + beta * (1 - sigma)) < lamp < 1,
        'c>0': c > 0,
        '0<delta<1/8': 0 < delta < F(1, 8),
        'eps(1+c)<1': eps * (1 + c) < 1,
        'eps+delta<1': eps + delta < 1,
    }
    # Gaussian line maps cost d p^(1/2+delta) alpha; paper alpha=Theta(p^(1/4+eps/2)),
    # package E alpha=ceil((8dp)^(1/4)) = Theta(p^(1/4+eps/4)).
    ga = F(1, 4) if pkg in ('E', 'F', 'G', 'H', 'I', 'J', 'K') else F(1, 2)
    cons['3/4+delta+(1+ga)eps<1'] = F(3, 4) + delta + (1 + ga) * eps < 1
    if pkg in ('A3', 'B'):
        cons['eps<1/12'] = eps < F(1, 12)
        cons['C1>=nu+3'] = C1 >= nu + 3
        cons['eps*C1<1'] = eps * C1 < 1
    else:
        if pkg in ('E', 'F', 'G', 'H', 'I', 'J', 'K'):
            # gamma = 2 d alpha^2 = O(b^(1/2+3eps/2)) = o(b); alpha < sqrt p needs eps < 1
            cons['eps<1/3 (gamma=o(b))'] = eps < F(1, 3)
        else:
            cons['eps<1/4 (gamma=o(b), alpha<sqrt p)'] = eps < F(1, 4)
        # guard: leaves have e<d^beta, at most (1-beta)log_m d+1 internal levels, <= m d pieces:
        # Delta = O(d^(2+beta+(1-beta)log_m B)) (crude piece count kept from the paper)
        cons['eps(2+beta+(1-beta)log_m B)<1'] = eps * (2 + beta + (1 - beta) * logmB_up) < 1
    if pkg == 'A3':
        cons['eps(2-tau)<1-tau'] = eps * (2 - tau) < 1 - tau
        g4 = 1 - tau - eps * (2 - tau)
    else:
        g4 = (1 - tau) * (1 - eps)
    g = dict(g1=1 - eps * (1 + c), g2=eps * c * (1 - tau), g3=eps * (1 - lamp), g4=g4,
             g5=F(1, 4) - delta - (1 + ga) * eps, g6=1 - delta - eps, g7=eps)
    bad = [k for k, v in cons.items() if not v]
    gmin = min(g.values())
    import math
    def l2(x): return math.log2(x)
    bh = BIT_STAGES if pkg in ('D', 'E') else (F_STAGES if pkg == 'F' else (G_STAGES if pkg == 'G' else (H_STAGES if pkg == 'H' else (I_STAGES if pkg == 'I' else (J_STAGES if pkg == 'J' else (K_STAGES if pkg == 'K' else bit_h))))))
    print(f"== package {pkg}: bit h={bh} complex h={cx_h}  nu={nu} log_m B<={float(logmB_up):.4f}")
    print(f"   1-tau = 2^{l2(th_b):.3f}  1-sigma = 2^{l2(th_c):.3f}  c = 2^{l2(c):.3f}  1-beta = 2^{l2(1-beta):.2f}  eps = {float(eps):.6f}  delta = {float(delta):.3g}")
    for k, v in g.items():
        print(f"   {k} = 2^{l2(v):.3f}")
    print("   all constraints hold" if not bad else f"   FAILED: {bad}")
    kappa = dyadic_below(gmin / 2, 8)
    assert 2 * kappa <= gmin
    print(f"   min g = 2^{l2(gmin):.3f};  kappa = min g/2 (paper convention) -> kappa = 2^{l2(kappa):.3f}; any kappa < 2^{l2(gmin):.3f} is valid")
    return dict(tau=tau, sigma=sigma, c=c, beta=beta, lam=lam, lamp=lamp, eps=eps, delta=delta, kappa=kappa, gmin=gmin, ok=not bad)


if __name__ == '__main__':
    out = {}
    for p in sys.argv[1:] or ('A3', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K'):
        out[p] = check(p)
