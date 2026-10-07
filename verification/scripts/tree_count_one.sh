#!/bin/bash
h=$1
python3 -c "h=$h; print(h); print(*range(h)); print(*[h-1-x for x in range(h)])" | "$(dirname "$0")/.."/csrc/tree > "$(dirname "$0")/.."/data/tree_$h.txt
