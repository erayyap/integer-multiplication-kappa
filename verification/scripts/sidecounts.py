import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import sys, json
from math import comb
sys.path.insert(0, f'{ROOT}/scripts')
from groupings import design, count
res = {}
for h in map(int, sys.argv[1].split(',')):
    for k in map(int, sys.argv[2].split(',')):
        gS, gT = design(h, k, "anti")
        tot, pairs, comp = count(h, gS, gT)
        v = comb(h, 3); z = 3 * comb(h - 3, 2)
        assert pairs == v * z
        res[f"{h},{k}"] = tot
        print(h, k, tot, round(tot / (v * z), 4), flush=True)
json.dump(res, open(f'{ROOT}/side_{sys.argv[1].replace(",","_")}.json', 'w'))
