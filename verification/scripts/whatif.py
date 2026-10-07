from math import comb, log
def th(h, lossf=1.0, zf=1.0, extra_w=0.0):
    v=comb(h,3); z=3*comb(h-3,2)*zf
    num=1-lossf*6*h*h/v
    if num<=0: return 0
    eta=num/((2+3*z+3*h/v+extra_w)*h**3)
    return -log(1-eta)/log(h**3)
for lossf in (1,0.5,0.25):
  for zf in (1,0.1,0.01,0.0):
    b=max((th(h,lossf,zf),h) for h in range(8,400))
    print(f"loss x{lossf:<5} side wires x{zf:<5}: best h={b[1]:3d} log2 theta={log(b[0],2):7.2f}  kappa~theta^2/6 = 2^{2*log(b[0],2)-2.585:.1f}")
