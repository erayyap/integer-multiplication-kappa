"""Optimize the fixed rational parameters of Section 8 for given motif exponents.

mode 'paper' : cost table exactly as in the paper (CRT/axis row d^2 l^tau).
mode 'linaxis': CRT/axis row d * l^tau (axes exposed by one non-adjacent swap each).
"""
import sys
from math import log2
from scipy.optimize import minimize
import numpy as np


def margins(eps, c, tau, sigma, beta, delta, mode, slack=1e-30):
    lam = max(tau, sigma, tau * (1 + c / beta))
    lamp = max(lam, sigma + beta * (1 - sigma))
    if lamp >= 1:
        return None
    g = dict(
        g1=1 - eps * (1 + c),
        g2=eps * c * (1 - tau),
        g3=eps * (1 - lamp),
        g4=(1 - tau - eps * (2 - tau)) if mode == "paper" else (1 - tau) * (1 - eps),
        g5=0.25 - delta - 1.5 * eps,
        g6=1 - delta - eps,
        g7=eps,
    )
    return g


def best(theta_t, theta_s, mode, C1=20):
    tau, sigma = 1 - theta_t, 1 - theta_s
    delta = 1 / 16
    epsmax = min(1 / 12, 1 / C1) * (1 - 1e-9)

    def f(x):
        le, lc, lb = x
        eps = min(2 ** le, epsmax)
        c = 2 ** lc
        beta = 1 / (1 + 2 ** (-lb))  # in (0,1)
        g = margins(eps, c, tau, sigma, beta, delta, mode)
        if g is None:
            return 1e9
        return -log2(max(min(g.values()), 1e-300))

    bestv, bestx = None, None
    for le0 in np.linspace(-60, -4.4, 12):
        for lc0 in np.linspace(-60, -2, 12):
            for lb0 in (-3, 0, 3, 10, 20, 30):
                r = minimize(f, [le0, lc0, lb0], method="Nelder-Mead",
                             options=dict(xatol=1e-6, fatol=1e-9, maxiter=4000))
                if bestv is None or r.fun < bestv:
                    bestv, bestx = r.fun, r.x
    le, lc, lb = bestx
    return -bestv, dict(eps=min(2 ** le, epsmax), c=2 ** lc, beta=1 / (1 + 2 ** (-lb)))


if __name__ == "__main__":
    tt, ts = float(sys.argv[1]), float(sys.argv[2])
    for mode in ("paper", "linaxis"):
        v, p = best(2 ** tt, 2 ** ts, mode)
        print(mode, "log2 min g =", round(v, 3), "-> kappa = 2^", round(v - 1, 3), p)
