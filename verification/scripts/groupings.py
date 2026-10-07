"""Generate source/target sunflower groupings and count side wires with csrc/blocks."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repo verification/ dir
import itertools, subprocess, sys
from math import comb

BIN = f"{ROOT}/csrc/blocks"


def grouping(h, k, prio, blocks):
    """prio: list rank per point (lower = preferred as petal).
    blocks: block id per point; groups = (kernel, block of petal), chunked to size<=k."""
    fam = {}
    for t in itertools.combinations(range(h), 3):
        p = min(t, key=lambda x: prio[x])
        ker = tuple(x for x in t if x != p)
        fam.setdefault((ker, blocks[p]), []).append(t)
    groups = []
    for key, ts in fam.items():
        ts.sort()
        for i in range(0, len(ts), k):
            groups.append(ts[i:i + k])
    return groups


def mask(t):
    return (1 << t[0]) | (1 << t[1]) | (1 << t[2])


def count(h, gS, gT):
    lines = [str(h), str(len(gS))]
    for g in gS:
        lines.append(" ".join([str(len(g))] + [str(mask(t)) for t in g]))
    lines.append(str(len(gT)))
    for g in gT:
        lines.append(" ".join([str(len(g))] + [str(mask(t)) for t in g]))
    out = subprocess.run([BIN], input="\n".join(lines), capture_output=True, text=True).stdout.split()
    count.extra = list(map(int,out)); return int(out[0]), int(out[1]), int(out[2])


def design(h, k, kind):
    if kind == "same":
        prio = list(range(h)); pS = pT = prio
    elif kind == "anti":
        pS = list(range(h)); pT = [h - 1 - x for x in range(h)]
    blocks = [x // k for x in range(h)]
    return grouping(h, k, pS, blocks), grouping(h, k, pT, blocks)


if __name__ == "__main__":
    h = int(sys.argv[1])
    for kind in ("same", "anti"):
        for k in map(int, sys.argv[2].split(",")):
            gS, gT = design(h, k, kind)
            tot, pairs, comp = count(h, gS, gT)
            v = comb(h, 3); z = 3 * comb(h - 3, 2)
            assert pairs == v * z
            print(f"h={h} {kind} k={k} wires={tot} ratio={tot/(v*z):.4f} z_eff={tot/v:.1f} complete={comp}", "cpairs,linew,linecells,sing=", count.extra[3:])
