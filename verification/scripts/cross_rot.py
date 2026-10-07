import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys
sys.path.insert(0,f'{ROOT}/scripts')
from cross_var import run, ORD
for h in map(int, sys.argv[1:]): print({'h':h,'rot_total':run(h, ORD['rot'])}, flush=True)
