import sys; sys.path.insert(0,'scripts')
from networks import counts, theta
from math import log
for h in list(range(18,40))+[46,100]:
    c=counts(h); W=c['Wc']; s=c['sc']; m=c['m']
    if c['eta_c']<=0: continue
    B=s+64*(W+m+1)**3
    nu=1
    while m**nu<B: nu+=1
    print(h, 'nu',nu,'C1=nu+3',nu+3, 'log2 theta_c %.2f'%log(theta(c['eta_c'],m),2), 'log2 theta_b %.2f'%(log(theta(c['eta_b'],m),2) if c['eta_b']>0 else float('nan')))
