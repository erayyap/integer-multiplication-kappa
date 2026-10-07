import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys, itertools
sys.path.insert(0,f'{ROOT}/scripts')
import fw_family as F
from math import log, comb
from fractions import Fraction as Fr
B={1:3,2:9.4}
res=[]
for k in range(2,9):
    for bits in itertools.product((0,1),repeat=k):
        fv=list(bits)+[1]; Z=[j for j in range(k) if fv[j]]
        if len(Z)!=1: continue
        j=Z[0]; r=k-j; beta=B.get(r,9.4)   # optimistic for r>=3
        sz=Fr(comb(k,j)*beta).limit_denominator(100)
        best=None
        for n in range(2*k,220):
            e=F.evaluate(n,k,fv)
            if not e: continue
            v,c,h=e['v'],e['c'],e['h']
            eta=Fr(v-6*c*h,(2*v+3*v*sz+3*c)*h**3)
            if eta<=0: continue
            th=-log(1-float(eta))/log(h**3)
            if best is None or th>best[0]: best=(th,n,c,h)
        if best: res.append((best[0],k,tuple(fv),best[1:]))
res.sort(reverse=True)
for th,k,fv,x in res[:10]: print(f'{log(th,2):.3f} k={k} f={fv} (n,c,h)={x}')
