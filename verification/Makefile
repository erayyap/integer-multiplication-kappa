CFLAGS ?= -O2
BINS = csrc/blocks csrc/gen csrc/grid csrc/sa_roles csrc/tree csrc/tree_bd csrc/tree_fake csrc/tree_rk

all: $(BINS)

csrc/%: csrc/%.c
	$(CC) $(CFLAGS) -o $@ $< -lm

verify: all
	python3 scripts/verify.py > data/verify_all.repro.txt
	diff data/verify_all.txt data/verify_all.repro.txt && echo "verify_all.txt reproduced exactly"

clean:
	rm -f $(BINS) data/verify_all.repro.txt
.PHONY: all verify clean
