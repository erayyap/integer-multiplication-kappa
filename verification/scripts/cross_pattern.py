import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys
sys.path.insert(0,f'{ROOT}/scripts')
import cross_width as cw, copy_sorted as cs
from collections import defaultdict, Counter
h=int(sys.argv[1]); n=h-1
st=cs.structure(n); mem=st['members']
srcn={g for g,_ in st['rects']}
ind=lambda t:[1 if i in t else 0 for i in range(h)]
TN={}
for q in range(h):
    o=[i for i in range(h) if i!=q]
    for nd in srcn: TN[(q,nd)]=frozenset(tuple(sorted((q,o[a],o[b]))) for a,b in mem[nd])
B={k:cw.rref([ind(t) for t in v],h) for k,v in TN.items()}
c={}
def leq(a,b):
    if (a,b) not in c: c[(a,b)]=TN[a]<=TN[b] or all(cw.inspan(ind(t),B[b]) for t in TN[a])
    return c[(a,b)]
by=defaultdict(set)
for k,v in TN.items():
    for t in v: by[t].add(k)
W={t:cw.dilworth(sorted(v),leq) for t,v in by.items()}
# per-kernel widths for comparison
pk={}
for t in by:
    pk[t]=sum(cw.dilworth(sorted(x for x in by[t] if x[0]==q),leq) for q in t)
print(Counter(W.values()), Counter(pk.values()))
for t in sorted(W):
    i,j,k=t
    if i in (0,1,2,h//2) and (k>=h-3 or k==j+1 or j==i+1 or 1):
        pass
# print matrix for i=h//3 interior: W as function of (j,k)
i=4
print('i=',i)
for j in range(i+1,h-1):
    print(j, ''.join(str(W[(i,j,k)]) for k in range(j+1,h)))
print('pk')
for j in range(i+1,h-1):
    print(j, ''.join(str(pk[(i,j,k)]%10) for k in range(j+1,h)))
