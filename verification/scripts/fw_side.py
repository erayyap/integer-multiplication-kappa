"""theta for FW families when side cost per element = (#kernels per element) * beta,
kernel = S cap T (unique when |Z|=1), beta = cover cost per element-kernel per (2 sides + rects)."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys
sys.path.insert(0,f'{ROOT}/scripts')
import fw_family as F
from math import log, comb
from fractions import Fraction as Fr
def best(k, fv, side_fn):
    out=None
    for n in range(2*k,200):
        r=F.evaluate(n,k,fv)
        if not r: continue
        v,c,h=r['v'],r['c'],r['h']
        sz=Fr(side_fn(n)).limit_denominator(1000)
        eta=Fr(v-6*c*h,(2*v+3*v*sz+3*c)*h**3)
        if eta<=0: continue
        th=-log(1-float(eta))/log(h**3)
        if out is None or th>out[0]: out=(th,n)
    return out
for beta in (3,5,9.4,12,15,25,40):
    a=best(3,[0,1,0,1],lambda n: 3*beta)
    b=best(7,[0,0,0,1,0,0,0,1],lambda n: 35*beta)
    c5=best(5,[0,0,1,0,0,1],lambda n: 10*beta) if True else None
    print(f'beta={beta:5}: k3 log2th={log(a[0],2):.3f}(n={a[1]})  k7 {log(b[0],2):.3f}(n={b[1]})')
