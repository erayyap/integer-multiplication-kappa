"""Numerical search: k-dim subspace orthogonal representations of the triple graph
(|S cap T| = 1  =>  T_S _|_ T_T) in R^f with an indefinite diagonal form of p positive signs.
Loss = sum_adj ||B_S^T E B_T||^2 + mu * sum_S ||B_S^T E B_S - G0||^2.
Prints final residuals.  usage: subrep.py h k f p [restarts] [g0sig]
"""
import sys, itertools, torch, math
torch.set_default_dtype(torch.float64)
h, k, f, p = map(int, sys.argv[1:5])
R = int(sys.argv[5]) if len(sys.argv) > 5 else 4
g0s = sys.argv[6] if len(sys.argv) > 6 else '+' * k
T = list(itertools.combinations(range(h), 3)); v = len(T)
idx = {t: i for i, t in enumerate(T)}
adj = [(i, j) for i in range(v) for j in range(i + 1, v) if len(set(T[i]) & set(T[j])) == 1]
I = torch.tensor([a for a, b in adj]); J = torch.tensor([b for a, b in adj])
E = torch.tensor([1.0] * p + [-1.0] * (f - p))
G0 = torch.diag(torch.tensor([1.0 if c == '+' else -1.0 for c in g0s]))
best = None
for r in range(R):
    torch.manual_seed(r)
    B = (torch.randn(v, f, k) / math.sqrt(f)).requires_grad_()
    opt = torch.optim.LBFGS([B], lr=1, max_iter=4000, history_size=50, tolerance_grad=1e-14, tolerance_change=1e-16, line_search_fn='strong_wolfe')
    def loss_fn():
        EB = B * E[None, :, None]
        cross = torch.einsum('afk,afl->akl', B[I], EB[J])
        self_ = torch.einsum('afk,afl->akl', B, EB) - G0
        return (cross ** 2).sum() + 1.0 * (self_ ** 2).sum()
    def closure():
        opt.zero_grad(); l = loss_fn(); l.backward(); return l
    for it in range(6):
        opt.step(closure)
    with torch.no_grad():
        EB = B * E[None, :, None]
        cross = torch.einsum('afk,afl->akl', B[I], EB[J]).abs().max().item()
        self_ = (torch.einsum('afk,afl->akl', B, EB) - G0).abs().max().item()
        l = loss_fn().item()
    print(f'restart {r}: loss {l:.3e} max|cross| {cross:.2e} max|self-G0| {self_:.2e}', flush=True)
    if best is None or l < best: best = l
print(f'h={h} k={k} f={f} p={p} ratio={f/k:.2f} best loss {best:.3e}')
