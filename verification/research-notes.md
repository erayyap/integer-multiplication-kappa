# Progress log (research phase, write-up after 07:30)

## Verified packages (scripts/verify.py, exact Fractions + interval logs)
- A3: kappa = 2^-109.08  (re-tuned constants, bit h=46, complex h=25)
- B : kappa = 2^-76.20   (O(d) non-adjacent axis swaps; C1 = 14)
- C : kappa = 2^-74.98   (delta tiny, eps<1/4 relaxation, beta-aware guard)
- D : kappa = 2^-73.73   (NEW 02:55) bit network: stage-1 scratch wires continued as
      stage-3 scratch wires with zero label loss; coordinate designs (48,42,48).
      theta_b = 2^-35.071 (was 2^-35.693).
- E : kappa = 2^-73.47   (NEW 03:05) alpha = ceil((8dp)^(1/4)) instead of
      ceil((12 d^2 b)^(1/4)); Gaussian row 3/4+delta+5eps/4 so eps < 1/5.
- F : kappa = 2^-70.14   (NEW 03:40) grouped side-wire gates (sunflower groups, size<=5,
      anti petal priorities); bit stages (48,44,48); theta_b = 2^-33.407.
- G : kappa = 2^-68.62   (NEW 03:50) laminar binary-tree groups per kernel family;
      stages (49,43,49); theta_b = 2^-32.644.  side/(v z) = 0.181 at h=48.

## F/G justification (side-wire compression)
Wire = source, gates touching it in time order, sink (paper l.534): extra gates on a
data wire add no wires.  Source (X) gates on the physical X wire of triple t:
own (B(x)F)+(P(x)<t>) < node_1 < ... < root (B(x)F)+(P(x)span t_root) < A(x)F.
Target (Y) gates: B(x)F < root' < ... < own'  with labels (B(x)F)+(P(x)(span t_node)^perp).
A side wire from source node G to target node H carries sum_{X in G} x_X and is added
to every Y_S, S in H; its chain BF < x(G) < y(H) < AF is nested iff span t_G _|_ span t_H
iff every (X,S) in GxH has |X cap S| = 1.  F2 requirement: centre + side = I, i.e.
side rectangles partition the neighbour pairs.  Node spans are sunflowers (Gram I+J).
Checked explicitly (scripts/check_invocation.py, scripts/check_tree.py, h=10).
12W gate bound of the paper is only used for the complex network (unchanged).
The orthonormal-basis requirement (paper l.789, 'binary label space') concerns the complex
network's F_2 labels only; the bit network uses F = Q^h with form I - J/9 where only
nondegeneracy is needed (sunflower spans: Gram I+J, positive definite).

## Side-wire experiments (04:00-04:10)
- tree DP (G): side/(vz) = 0.181 (h=48), 0.152 (h=64).
- breakdown h=40: KK wires 0.26M (7.2M cells), disjoint-kernel 3.73M wires (12.5M cells),
  of which PP singletons 0.66M; defects (|X cap S|=2 inside KP/PK columns) cost ~1.4M.
- 'same' priorities: KP/PK defect-free but PP = 1/5 of pairs -> worse (5.88M vs 3.99M).
- rank-based petals (min/mid/max x min/mid/max): anti (min,max) best.
- grid / slice hierarchies merging families (csrc/grid.c, csrc/gen.c): best S1A x S2B,
  only 1.1-1.7% fewer wires.  SA over per-triple roles (csrc/sa_roles.c): ~1-2%.
- Within-stage scratch reuse impossible: invocation scratch label is forced up to
  F^{(x)j} (x) Q; only stage1 -> stage3 nests (D).  Contaminating gates costs ~h loss per data wire.
- eps < 1/5 is forced by eta = 1/(4d) prime-interval requirement (theta_i ~ 1/d).

## D justification (scratch reuse 1 -> 3)
Stage-1 wire (invocation a2,a3) ends with label F (x) t_{a2} (x) t_{a3}
(time 7 / time 6 gate label (A(x)F)(x)Q, A trivial).
Stage-3 invocation (a1',a2'): A=F(x)F, P=t_{a1'}(x)t_{a2'}, B=P^perp, Q trivial.
 side wire first gate (time 0, target Y with triple S'): label B (x) <t_{S'}>.
 center first gate (time 1): B (x) F.
If |a2 cap a2'| = 1 (t_{a2} perp t_{a2'}) then F(x)t_{a2} is inside B, so
 F(x)t_{a2}(x)t_{a3} is inside B(x)<t_{S'}> when S'=a3, and inside B(x)F always.
Continuation edge increases dimension; residual = orth. complement, nondegenerate.
Matching: z2-regular bipartite group graph, equal supplies => integral b-matching.
Stage 2 cannot be chained (loss 1 per wire from the P-component; ordering
argument shows X-first side pattern creates a cycle with the center order).

## E justification
Lemma 7.1 needs theta_i = t_i/s_i - 1 > p/alpha^4 and 2 <= alpha < sqrt p.
theta_i > eta = 1/(4d).  alpha^4 >= 8dp gives alpha^4 theta_i > 2p.
gamma = 2 d alpha^2 = O(d^{3/2} p^{1/2}) = o(b) for eps < 1/3 (final error needs gamma <= b/4).
Cost d p^{1/2+delta} alpha = p^{3/4+delta+5eps/4}.

## Dead ends this phase
- F_q labels for the bit network are allowed (digitwise radix-q frames) but the
  triple graph |S cap T|=1 contains Kneser K(n-1,2); Haviv's topological bound
  gives orthogonality dimension >= n-3 over every field: no gain.
- Asymmetric coordinates without reuse: symmetric optimal.
- General bound: eta <~ (1-6/R)/(3 z xi^3) with z >= R xi/2 ; triples are within
  constant factors of what such designs can do at their R.

## Package H (copy wires, 04:20) : kappa = 2^-61.114
Per kernel point q the neighbour cells (q+P, q+Q), P cap Q = {} (pairs of the other n=h-1 points)
are partitioned by rectangles G x H (supports disjoint) = side wires x(G)=BF+P(x)span t_G ->
y(H)=BF+P(x)(span t_H)^perp.  Source copies of x_P (copy wires) and target accumulators carry
chains of nested nodes; #chains per element = Dilworth width of its nodes.  Data wire X (Y) is
the first source (target) chain of family q=min.  side(h) = h*total(h-1) - 2 C(h,3).
Garbage handling: s -= sum a at BF before the copy; restores at AF (see check_copy.py).
scripts/check_copy.py builds a full stage-1 invocation (all wires, gates, labels over Q^h mod p)
and checks F2 map y += x with arbitrary scratch restored, nesting of every label chain
(only decreasing edges: centre 3->4), nondegeneracy of all gate labels.  PASS h=8,10.
Tree rule (copy_cover.py): stages (51,48,51), 1-tau = 2^-28.894.

## Package I (sorted cover, 04:40) : kappa = 2^-60.017
scripts/copy_sorted.py: rule by the sorted pattern of P=(x1<x2), Q=(z1<z2) (PPQQ, QQPP, PQPQ,
QPQP, PQQP, QPPQ) with nodes Wpre/rowfrom/rowsuf/colpre/colmid; ~3.74 C(n,2) rects, width 3.
total(47) = 10166 vs 14594 (tree).  check_copy PASS h=8,10 (copy_sorted).  Stages (51,50,51),
1-tau = 2^-28.345, g2 binding.

## Structural facts (04:50)
- Three stages are forced: the endpoint identity needs rho(X_a)=Y_a (a swap) = 3 shears, each
  stage exposing one tensor factor; so m = h1 h2 h3.  Prop 11 works for any m with s<Wm.
- Centre: every centre wire loses exactly h (P(x)F) whatever it gathers, so L/I = c*h with
  c = min F2-rank of C (diag 1, zero on |S cap T| in {0,2}, free on 1).  Pencil {a,b}: C restricted
  = I_{h-2}, so c >= h-2 (paper: c=h).  Point-linear factorizations force c=h.
- Sensitivity (theta_b, log2): side wires x0.5 -> +0.55 (kappa +1.1 per halving... ~+1.9 at x0.5
  incl. re-optimized h); centre loss x0.5 -> theta +3.0 (h drops to ~26-31).
- Centre min-rank (CP-SAT, center/minrank.py): h=6 c=4 infeasible; h=7 c=6 feasible (uses two
  Fano planes, does not scale); h=8 c=6 infeasible.  So c >= h-1 in practice; gain <= 1/h of L.
- Refined centre labels (gather label = star hyperplane) save at most ~2 per invocation (nesting
  along X/Y forces the other centres to full labels).  Negligible.
- Relaxed side labels: node label U(A) for any point set A with supp G <= A, A cap B_H = {};
  chains only need nested supports; F2 parity => odd covers suffice.  check_copy.py supports
  'gsupp'/'hsupp' (sorted_supp.py PASS h=8, bad_supp.py FAILS as it should).
  For copy_sorted the support widths equal the member widths (its label intervals are tight).
- F2 rank of the pair-disjointness matrix (relax/rank_lb.py): = C(n,2) for n=8,12,20 (44/45 at 10);
  drops for n = 9, 11.  So rects per kernel >= ~C(n,2) even with odd covers.
- check_copy copy_sorted PASS at h=11 (scratch 3245) and h=12 (4492), = formula h*total(h-1)-2v+h.
- Stage 2 = dual network: gates reversed (F2 XOR gates are involutions), banks renamed X<->Y,
  every local label U -> U^perp (inclusions reverse, time reverses => chains nested again; centre
  loss preserved; endpoints map X:<t>->F to Y:0->t^perp and back).  check_copy.py DUAL=1 builds it
  explicitly: copy_sorted h=8 PASS (x += y, all restored, nested, nondegenerate).
  Stage j>=2 embedding: local U -> (B(x)F) + (P(x)U) is monotone and nondegenerate; all scratch
  first labels are BF, last AF, so the 1->3 reuse (D) is unchanged.
- DUAL h=10: copy_sorted PASS (2490 wires nested, scratch 2250) and copy_cover PASS (2790 wires, scratch 2550);
  both give x += y with every scratch restored.

## 05:30-05:45 explored and discarded
- Transform-layout recursion (wide chunk swaps instead of D*q width-K swaps): a breakpoint/run
  argument (each swap fixes <= 4 boundaries; runs start at width K) shows ~Dq/8 swaps of width >= K
  are needed, so the layout cost p K^{tau-1} is optimal for contiguous-block swaps.  The layer's
  real condition is lambda' > tau(1+c) (not tau(1+c/beta)), but beta ~ 1 already, so no gain:
  kappa ~ eps*theta^2/2 is structural (layout K^{-theta} vs layer (dK)^tau).
- Block centres (G = H = triples inside Y, loss |Y|): nesting along X_S forces the blocks
  containing S to be a chain B1 < B2 < ...; 4-set parity (|S cap T| = 2 pairs even) forces every
  point of B1 \ S to make an odd suffix, so B1 = S: every triple becomes its own centre (loss v).
  Dead.  Complement-star blocks (loss h(h-1)) are incomparable on X_S, so also invalid.
- Interface saving: endpoint rank of a data pair is rank(I+D)+rank(I-D) >= m with D = M_in(Y)-M_in(X);
  the paper (D = P_U) pays 2m-1, i.e. saves only 1 per pair.  A pair with Y source -P_V (V perp U)
  and X sink P_{V^perp} would save dim V, but needs V perp every X label (X reaches A(x)F at the
  stage-3 centres) and V inside Y's first label (in F(x)t2(x)t3, which X already covers): only
  isotropic V could work, and projections need nondegenerate labels.  No construction found.
- Constants: theta_c (complex) is far from binding (needs theta_c < theta_b^2); eps < 1/5 from the
  Gaussian row (only +0.74 bits even if relaxed to the gamma bound 1/3).  Only theta_b matters.
- Point-indexed centres G_i, H_i inside star_i with order-type rules g_i(S) = g[rank of i in S],
  h likewise: the |S cap T| = 2 conditions for the three positions of the third point give
  g_min h_min = g_mid h_mid = g_max h_max, and the diagonal forces the common value 1, so
  g = h = (1,1,1), which is the paper's star centre.  Min-stars (each triple once, loss ~h^2/2, nested)
  violate the |S cap T| = 2 parity.  Pricing: 1 unit of loss ~ 360 wires at package I, so centre loss
  (h^2 per invocation) and side wires (26 v) carry comparable weight.
- Cycle bound on terminals: for every rho-cycle, sum of endpoint ranks >= m (rank of k*I), so in principle
  a k-cycle of scratch roles can cost ~m instead of k*m.  With projection labels, wire w with source
  -P_U and predecessor ending at U^perp saves dim U; the paper uses this only for X_a (dim U_a = 1).
  For a scratch wire p ending at L_p (stage-2 wires end in F(x)F(x)t_{a3}) and successor w starting at
  G_w, saving = dim(G_w cap L_p^perp); but realizing rho(p) = w needs value routing through shared
  gates whose labels would push p up to F again.  Not realized; open.

## 05:30 restore-at-FULL; side-cover fork result; headroom
- Restore gates (time 7.x) read X after centre gather (time 3) => scratch ends at A(x)F label. Continuing stage1->2 or 2->3 costs loss 1/wire (~1/360 of a unit wire at pkg I) -> worse. Only stage 1->3 continuation (package D) is free.
- Fork agent (cover/): no cover beats copy_sorted. Per-order-type F2 rank = C(n-2,2) for each of the 4 interleaved types => sorted (4048 rects at n=47) is within 2.2% of the per-type bound 3960. Small-n CP-SAT: n=5 40 vs 44, n=6 76 vs 80 (boundary effects). Product/tree covers far worse. |G|=1 direct-scatter trick blocked by chain order.
  Open: joint rank of interleaved types ~C(n,2) (4x smaller) but needs node families beyond prefix/suffix/interval supports.
- Headroom (scripts/headroom.py), stages (51,50,51) optimal throughout:
  side/v 26.2 -> kappa 2^-60.00 ; 13.1 -> 2^-58.11 ; 7.0 -> 2^-56.47 ; 3.0 -> 2^-54.47 ; 0 -> 2^-50.48.
  So even free side wires give at most +9.5 bits over package I; realistic cover lower bound (~7) gives +3.5 bits.
- 05:35 check_copy copy_sorted h=13 (stage 1) PASS; DUAL h=12 (stage 2) PASS. side+copies = 1920+2560 = 4480 (h=12), 2574+3432 = 6006 (h=13) == h*total(h-1) - 2C(h,3) formula exactly.
- 05:50 FW families with copy-cover sides (scripts/fw_side.py, fw_side2.py): side/elt = C(k,j)*beta_{k-j}
  (kernel = S cap T unique when |Z|=1). Only k=7, f=C(j,3) mod 2 (Z={3}, centres = triples, n=28) competes:
  beats triples (beta_2=9.4 achieved) iff beta_4 <= ~12.4 for 4-set disjointness covers on 25 pts.
  Sorted-pattern generalisation to 4-sets needs ~70 patterns x C(n,6)/C(n,4) ~ 1000 rects/elt -> dead.
  k=6 Z={2}: -29.69 even at beta 9.4 (worse). No |Z|>=2 family relevant.
- 05:38 copy_sorted h=14 stage-1 PASS (8582 wires nested, 4282 gate labels nondegenerate), DUAL h=13 stage-2 PASS.
  verify.py re-run reproduces data/verify_all.txt byte-for-byte.

## Package J (05:50) : kappa = 2^-59.591  -- cross-kernel copy chains
Loophole in the old "cross-kernel chains impossible" claim: that was about side-wire orthogonality
(<t_S,t_T> = -1 for disjoint triples).  Copy wires / accumulators only need NESTED spans along
the wire.  A kernel-q node that is a star at a (all triples contain {q,a}) has span inside a
kernel-a node containing the same triples (equal spans allowed: label stays constant).
scripts/cross_width.py: per triple, Dilworth width of its source nodes over all 3 kernels under
span inclusion (exact mod p), global tie-break by (span rank, key) => one global time order.
Interior triples: width 7 per side (was 3+3+3 = 9).  cross(h) = (7h^3 - 39h^2 + 26h + 102)/6
exact for h = 9..20 (fit on 14..20, predicts 9..13,15 exactly).  side = h*rects(h-1) + 2 cross(h) - 2v.
scripts/check_cross.py (global chains, gates sorted by span rank): h=8 PASS stage 1 and DUAL.
Stages (51,50,51) theta_b = 2^-28.130 (I: -28.345).  verify.py J -> kappa = 2^-59.591.
- check_cross (global order) PASS h=8,10,12 stage 1 + DUAL; counts = formula.

## Package K (06:00) : kappa = 2^-59.335 -- cyclic per-kernel orders
Kernel q orders the other points q+1, q+2, ..., q-1 (mod h) (scripts/cross_var.py 'rot').
Total cross-kernel chains (src+tgt) = 2h(h-3)^2 EXACT for h = 9..19 (data/cross_rot_*.log):
interior width 6 per side (global order 7, per-kernel 9).  Rule variants (PQQP/QPPQ A/B) identical.
check_cross.py ORDER=rot: PASS h=8,10,12 (stage 1 and DUAL); copies = 2h(h-3)^2 - 2v exactly.
(h=9 fails nondegeneracy: I - J/9 itself is degenerate at h=9; irrelevant.)
Stages (51,50,51), theta_b = 2^-28.004, verify.py K -> kappa = 2^-59.335.
Other orders tried (h=12/13): reverse, altrev, belowRevFirst, rotmid, rotalt, zigzag: worse.
- 05:50 rect merging (scripts/merge_rects.py): rects sharing a source (target) triple set -> one wire
  with union target (source) node; greedy accepting only total decreases: h=10 1940 -> 1900 (-2%).
  Union nodes create antichains in the other side's chains, so most merges don't pay.
- Chain anatomy under rot (h=16, triple (3,8,12)): per side 3 colpre/Wpre chains (one per kernel)
  + 3 chains each merging kernel-x rowsuf with kernel-y colmid stars.  colpre stars (leaves q+1,q+2,..)
  would need 'rowmid' stars in the centre's kernel; rule variants providing them give the same total.
- 05:52 K confirmations: rot formula 2h(h-3)^2 also exact at h=23 (18400); global J cubic exact at h=26 (16241).
  check_cross ORDER=rot h=14 PASS stage 1 + DUAL (6762 wires nested, 4006 labels nondegenerate, copies 2660 = formula).
  Variants under rot: nested-type A/B (4 combos), separated-type keys (z1 vs x2; x1 vs z2): equal or worse;
  rotrev = rot; rotalt worse.  verify_all.txt regenerated (A3..I unchanged; J -59.591; K -59.335).
- 06:05 twin rects: under rot + QPPQ-B, kernel-q PQPQ rect (colpre(d,c) x rowsuf(d,c)) and kernel-c QPPQ-B rect
  (rowpre(q,d) x colpre(q,d)) have the SAME source triple set {q,c,y: y in (q,d)} and the same target leaf set
  (c,q); one side wire G x (H u H') is valid (|S cap T| = 1 throughout).  But the union target nodes are
  incomparable with the targets' other nodes: h=10 merging all 150 twins: rects -150, chains +590 (net worse).
  Greedy single merges: AA/AB/BA -20..-40 only (2%).  Dead end unless targets are redesigned around unions.
- J cubic also exact at h=30 (25797), h=34 (38505); K formula exact at h=26 (27508).
- 06:10 rot width anatomy (scripts/rot_gaps.py, h=14, src and tgt identical): triple with cyclic gaps (g1,g2,g3):
  all >= 2 -> 6; exactly one gap = 1 -> 3 (2 for exactly one composition per position);
  two gaps = 1 -> 2.  Sum = (h/3)[6 C(h-4,2) + 3(3(h-5)+2) + 6] = h(h-3)^2 per side  (matches the fit).
- 06:12 per-kernel mixing of rule variants under rot (scripts/cross_mix.py, 6 schemes): all 2794 at h=11 (= K). 6/side appears structural for these node families.
- 06:08 check_cross ORDER=rot h=16 PASS stage 1 + DUAL (10800 wires nested, 6402 labels nondegenerate, copies 4288 = formula).
- K formula exact at h=31 (48608).
- 06:33 SA over independent per-kernel orders (scripts/sa_orders.py, 25 min x 3 at h=10, 1 at h=8; ~17k + 30k evals): no improvement on rot (980 / 400). Cyclic-offset hill-climb (offset_search.py h=10,12): none so far. rot appears optimal for this cover.
- 07:06 cyclic-offset hill-climb finished (80 min each, h=10 and h=12): final = rot (980, 1944); no better cyclic order.
