# Improving κ in integer multiplication below n log n

This repository contains a manuscript and verification code that improve the exponent κ in the O(n (log n)^{1−κ}) bound of the OpenAI preprint *Integer multiplication below n log n*, from 2^-182 to 2^-59.335. Its layout follows [openai/math](https://github.com/openai/math).

The results were produced by **Claude Opus 5.5** (Anthropic), working as an agent in Claude Code sessions directed by the repository owner. They are computer-assisted and have **not been refereed**. None are formalized in Lean. Some of the arguments are given in outline only; each manuscript says exactly which claims are machine-checked and which are not. Corrections are welcome via issues.

## Navigating the collection

- Use the [manuscript map](CONTENTS.md) to find individual papers and their supporting materials.
- The [`preprints/`](preprints/) directory contains PDFs, LaTeX sources, and per-manuscript citation and build instructions.
- The [`verification/`](verification/README.md) directory contains the exact verifiers, gate-level network checkers, search code, logs and research notes.

## Current contents

| # | Result | Paper |
|---|---|---|
| 001 | Improves the exponent in OpenAI's *Integer multiplication below n log n* from κ = 2^-182 to **κ = 2^-59.335** | [PDF](preprints/Improving-kappa-in-integer-multiplication-below-n-log-n-October-7-2026/paper.pdf) |

## How the results were produced

Result 001 started from OpenAI's preprint and asked for a larger κ. The model read the source paper, wrote exact verifiers for its constraint system, searched for better constructions (combinatorial covers, label designs, parameter changes), and checked every numerical claim with exact rational arithmetic. New networks were checked gate by gate for small parameters. The written manuscript was generated from the research notes, which are kept unedited in [`verification/research-notes.md`](verification/research-notes.md).

## Versions and citations

Corrections and revisions will be recorded as new commits, with earlier versions kept in the git history. To cite an individual manuscript, use the BibTeX block in its directory.

## License

Apache License 2.0, see [LICENSE](LICENSE).
