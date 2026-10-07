import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys
sys.path.insert(0,f'{ROOT}/scripts')
import cross_width as cw, copy_sorted as cs
from collections import defaultdict, Counter
h=int(sys.argv[1]); n=h-1
st=cs.structure(n); mem=st['members']
for role,nodes in (('src',{g for g,_ in st['rects']}),('tgt',{x for _,x in st['rects']})):
    ind=lambda t:[1 if i in t else 0 for i in range(h)]
    TN={}
    for q in range(h):
        o=[(q+1+i)%h for i in range(h-1)]
        for nd in nodes: TN[(q,nd)]=frozenset(tuple(sorted((q,o[a],o[b]))) for a,b in mem[nd])
    B={k:cw.rref([ind(t) for t in v],h) for k,v in TN.items()}
    c={}
    def leq(a,b):
        if (a,b) not in c: c[(a,b)]=TN[a]<=TN[b] or all(cw.inspan(ind(t),B[b]) for t in TN[a])
        return c[(a,b)]
    by=defaultdict(set)
    for k,v in TN.items():
        for t in v: by[t].add(k)
    tab=defaultdict(set)
    for t,v in by.items():
        if 0 not in t: continue
        a,b,cc=t; g=(b-a,cc-b,h-cc+a)
        w=cw.dilworth(sorted(v,key=lambda k:(len(B[k]),str(k))),leq)
        tab[g].add(w)
    print(role)
    cnt=Counter()
    for g in sorted(tab):
        cnt[(tuple(x==1 for x in g), tuple(sorted(tab[g])))]+=1
    for k,v in sorted(cnt.items()): print('  gap==1 pattern',k[0],'widths',k[1],'count',v)
