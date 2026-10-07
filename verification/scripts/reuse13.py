# Stage-1 scratch wires reused (zero loss) as stage-3 scratch wires.
# Stage-1 final label F (x) t_{a2} (x) t_{a3} lies inside stage-3 side label
# B(x)<t_Y> (B=(t_{a1'}(x)t_{a2'})^perp) when t_Y=t_{a3}, t_{a2} perp t_{a2'}.
from math import comb, log
from fractions import Fraction as Fr
def eta(hs, reuse=True):
    N=1; m=1
    for h in hs: N*=comb(h,3); m*=h
    per=[]; Ln=Fr(0)
    for h in hs:
        v=comb(h,3); z=3*comb(h-3,2); c=h
        per.append(z+Fr(c,v)); Ln+=Fr(c*h,v)
    Wn = 2 + (max(per[0],per[2])+per[1] if reuse else sum(per))
    return (1-2*Ln)/(Wn*m), m
for reuse in (False,True):
    best=max(((-log(1-float(eta(hs,reuse)[0]))/log(eta(hs,reuse)[1]),hs)
              for a in range(30,70) for b in range(30,90) for hs in [(a,b,a)] if eta(hs,reuse)[0]>0))
    print("reuse" if reuse else "paper", best[1], "log2 theta_b =", round(log(best[0],2),3))
