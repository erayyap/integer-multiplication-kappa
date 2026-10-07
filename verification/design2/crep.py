"""Complex symmetric minrank search for the triple graph (|S cap T| = 1 => t_S^T t_T = 0),
t_S in C^f, t_S^T t_S = 1 (complex bilinear, not Hermitian).  A solution over an algebraic
number field K of degree d gives a d-dim subspace rep over Q in Q^{d f} (xi = f).
usage: crep.py h f [restarts] [kind]   kind: tri (exactly-1 graph)"""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys, itertools, torch, math
torch.set_default_dtype(torch.float64)
h, f = int(sys.argv[1]), int(sys.argv[2]); R = int(sys.argv[3]) if len(sys.argv) > 3 else 4
T = list(itertools.combinations(range(h), 3)); v = len(T)
adj = [(i, j) for i in range(v) for j in range(i + 1, v) if len(set(T[i]) & set(T[j])) == 1]
I = torch.tensor([a for a, b in adj]); J = torch.tensor([b for a, b in adj])
best = None
for r in range(R):
    torch.manual_seed(r)
    X = (torch.randn(v, f) / math.sqrt(f)).requires_grad_(); Y = (torch.randn(v, f) / math.sqrt(f)).requires_grad_()
    opt = torch.optim.LBFGS([X, Y], lr=1, max_iter=5000, history_size=50, tolerance_grad=1e-15, tolerance_change=1e-18, line_search_fn='strong_wolfe')
    def loss_fn():
        re = (X[I] * X[J] - Y[I] * Y[J]).sum(1); im = (X[I] * Y[J] + Y[I] * X[J]).sum(1)
        sre = (X * X - Y * Y).sum(1) - 1; sim = 2 * (X * Y).sum(1)
        return (re ** 2 + im ** 2).sum() + (sre ** 2 + sim ** 2).sum()
    def closure():
        opt.zero_grad(); l = loss_fn(); l.backward(); return l
    for it in range(8): opt.step(closure)
    l = loss_fn().item()
    print(f'restart {r}: loss {l:.3e}', flush=True)
    best = l if best is None else min(best, l)
    if l < 1e-13: torch.save((X.detach(), Y.detach()), f'{ROOT}/data/subrep/crep_h{h}_f{f}_r{r}.pt')
print(f'h={h} f={f} best {best:.3e}')
