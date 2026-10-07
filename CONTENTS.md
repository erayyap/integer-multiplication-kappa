# Manuscript map

<table>
<thead><tr><th>Result</th></tr></thead>
<tbody><tr>
<td>

**001. Improving κ in integer multiplication below n log n.** OpenAI's preprint gives a fixed multitape Turing machine multiplying n-bit integers in time $`O(n(\lg n)^{1-\kappa})`$ with $`\kappa = 2^{-182}`$. Keeping its architecture, this result raises the exponent to $`\kappa = 2^{-59.335}`$. The gain comes mainly from a new bit-level motif network (rectangle side wires fed by copy chains shared across kernels, with cyclic point orders) and from re-tuned constants. Computer-assisted; not refereed.

</td>
</tr></tbody>
<tbody><tr>
<td>

&emsp;[Integer multiplication below n log n: improving κ from 2^-182 to 2^-59.335](preprints/Improving-kappa-in-integer-multiplication-below-n-log-n-October-7-2026/paper.pdf)

The OpenAI preprint *Integer multiplication below n log n* gives a deterministic multitape Turing machine that multiplies two n-bit integers in time $`O(n(\lg n)^{1-\kappa})`$ with $`\kappa = 2^{-182}`$. We keep the architecture of that algorithm and improve the exponent to $`\kappa = 2^{-59.335}`$. Every package is checked by an exact rational verifier, and the new bit network is built gate by gate and checked for $`h \le 16`$. We also record structural limits (this architecture cannot pass about $`2^{-50.5}`$) and negative results on alternative designs.

Supporting material: [`verification/`](verification/README.md)

</td>
</tr></tbody>
</table>
