from math import comb, log
def stage(n): return comb(n,3), n, n, 3*comb(n-3,2)   # v,c,h,z
def eta(ns, reuse):
    S=[stage(n) for n in ns]; m=1
    for v,c,h,z in S: m*=h
    deficit = 1 - 2*sum(c*h/v for v,c,h,z in S)
    if deficit<=0: return 0
    w=2.0
    for j,(v,c,h,z) in enumerate(S):
        later=1
        for (_,_,hh,_) in S[j+1:]: later*=hh
        f = 1/later if reuse else 1
        w += (z + c/v)*f
    return deficit/(w*m), m
best={}
for reuse in (False,True):
    b=(0,)
    for n1 in range(30,90):
      for n2 in range(30,90):
        for n3 in range(30,90):
          r=eta((n1,n2,n3),reuse)
          if r and r[0]>0:
            th=-log(1-r[0])/log(r[1])
            if th>b[0]: b=(th,(n1,n2,n3))
    print('reuse' if reuse else 'noreuse', b, log(b[0],2))
