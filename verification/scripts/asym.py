# Asymmetric stages: coordinate j uses the triple design on h_j points.
from math import comb, log
from fractions import Fraction as Fr
def eta(hs, cplx=False):
    N=1
    for h in hs: N*=comb(h,3)
    m=1
    for h in hs: m*=h**3//h**2  # m = prod h_j
    m=1
    for h in hs: m*=h
    Wn=Fr(2); Ln=Fr(0)
    for h in hs:
        v=comb(h,3); z=3*comb(h-3,2) if not cplx else comb(h-3,3)+3*(h-3); c=h if not cplx else h+1
        Wn+=z+Fr(c,v); Ln+=Fr(c*h,v)
    num=1-2*Ln - (1 if cplx else 0)*0
    return num/(Wn*m), m
best=[]
for a in range(10,80):
  for b in range(a,90):
    for c in range(b,120):
      e,m=eta((a,b,c))
      if e<=0: continue
      th=-log(1-float(e))/log(m)
      best.append((th,(a,b,c)))
best.sort(reverse=True)
for th,hs in best[:8]: print(hs, "log2 theta_b =", round(log(th,2),3))
e,m=eta((46,46,46)); print("sym 46:", log(-log(1-float(e))/log(m),2))
