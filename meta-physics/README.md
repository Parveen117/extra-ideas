# Meta physics: certified developments from the idea vault

Research owner: **Monty Dabas**. Started 3 October 2026.

This folder develops ideas from Universal-Assurance-Ledger into scoped,
proved and reproducible results. The private PDFs stay at their source.
Only newly written developments, their source audit, and certificates
are published here. Yang–Mills development is paused for this work.

**Reading completed:** 98 PDFs, 299 pages, 95 distinct extracted texts.
The [full catalogue and queue](SOURCE_AUDIT.md) accounts for every file;
the [inventory](SOURCE_INVENTORY.json) retains source hashes and status.
Reading completion is not certification of all source claims.

## First certified development: MP-1

[Hidden frame memory and a certified observation protocol](mp1/THEOREM.md)
develops the winding-memory and observation-axis ideas, using the existing
native EMK carrier. Seven written results establish:

- The quadratic axis observer identifies g and -g. One complete observed
  turn reverses the hidden frame; two return it.
- An exact rational winding algorithm rejects zero-crossing and open paths.
- A computable distance-to-zero margin protects the count against every
  perturbation within a declared error tube, including between-sample error.
- Frame recovery through this observer needs a strict quarter-turn
  increment bound. At the boundary, two opposite frame answers fit the
  same observation.
- One parity record recovers the frame sign. Full winding needs a larger
  integer ledger. A retained reference can reveal the hidden sign.
- Global sheet memory can coexist with locally flat transport; it is not
  automatically non-Abelian curvature or entropy production.

These mechanisms have standard mathematical antecedents and existing
native predecessors. The proof gives the attribution and the precise
extension beyond RKF T76 and Extra Ideas R7/R14/R15/R39/R42.

The concrete one-turn diamond has winding 1 and squared clearance 1/2.
Vertex error 1/8 plus chord error 1/8 in each coordinate gives a squared
tube bound of 1/8, safely below 1/2. A boundary perturbation reaching
zero is explicitly rejected. No experimental noise distribution is assumed.

## What certification means here

Each result needs explicit hypotheses, a written proof, exact executable
controls and a source/dependency record. The [MP-1 certificate](mp1/VERIFICATION.json)
replays **23 tests**, including **6,561** exact vertex perturbations and
failure controls. Its general tolerance claim rests on the proof, not
on the number of examples. The workflow uses **Python 3.12 only** and
the standard library. It imports the existing rational EMK implementation
without altering the canonical native engine.

This is mathematical/code certification within the stated model. It is
not proof-assistant verification, external peer review or empirical
validation. Source 017/078 are partially developed; their full proposed
physical interpretations remain open. R1–R46 retain their original
register and evidence; MP numbers are a separate development sequence.

Replay from the repository root:

```sh
python3.12 -B meta-physics/mp1/verify.py --check
```

This checks the inventory, preserved R1–R46 register, dependency hashes,
all MP-1 tests and exact agreement with the saved certificate. No private
repository access is required. For verbose individual test output:

```sh
python3.12 -B -m unittest discover -s meta-physics/mp1 -p 'test_*.py' -v
```

## Next certification gates

1. Supply an actual native source/observer protocol whose between-sample
   error budget can be derived, completing MP-1's measurement contract.
2. Repair the NRT loop-observation family with identifiability probes and
   a distinction between signed circulation and nonnegative dissipation.
3. Develop the remaining PDF families through shared native results;
   duplicate drafts reference the same certified chapter.

No raw draft is promoted merely by changing its notation or moving it
into this folder.
