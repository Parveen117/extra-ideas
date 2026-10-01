# EMK constructive cut foundation

This standalone corrective edition develops **14 proved mathematical statements from one explicit typed cut-record construction**, with a source-pinned computational certificate. The original manuscript's **48 claims are all accounted for**. The original manuscript is **not certified as written**: its universal raw closure claim has an exact counterexample. A corrected construction does not retroactively prove that claim.

Read [the manuscript](emk_topology_foundation.pdf), [its TeX source](emk_topology_foundation.tex), and [the certificate](certificate/VERIFICATION.json). Extra-work applications and physical identifications are deferred. This folder has no dependency on earlier extra-ideas runtime code.

The construction proceeds from typed cut/self-cut records to signed counts and refinement, then to a free ledger on the two cut roles. Role exchange `K` and role parity `H` are constructed on that ledger. Their composite `R = KH` satisfies `R² = −I` and `R⁴ = I` by its action on the two free generators. **Neither an ordinary complex field nor a four-element phase alphabet is an input.** The full noncommuting operator algebra is retained. Its commuting scalar sector, completed only afterward, contains the derived iota. It does not replace full UGD states or their seam memory.

This is an existence/construction result, not a uniqueness theorem for arbitrary acts called cut and self-cut. Assigning the free tags `k → H`, `c → K` is an explicit representation: `ck` acts as `R`; four such composites contain eight primitive tags. The upload's claimed four-primitive-tag identity is not restored.

| Statements | Result | Evidence |
|---|---|---|
| C1–C3 | Typed records, countermodel to raw closure, signed cut refinements | Written proofs; exact countermodel |
| C4–C7 | Derived quarter-turn, full operator algebra, scalar sector, pairing, seam aperture, mixed commutator and cut-corner residue | Written proofs; complete generator/bilinear basis checks; canonical RKF replay |
| C8 | Complete coefficient field, radial order, native square root and norm | Written Cauchy and bisection proofs |
| C9–C12 | Recognition topology, prefix ultrametric, complete finite/infinite history space and correctly typed Eye hierarchy | Written general proofs; finite exhaustive checks |
| C13–C14 | History retained at operator return; exact path-equality criterion including global loop periods | Written proofs; exact graph checks |
| N1–N4 | Euler exponential, positive logarithm, native period and scalar polar/logarithm | Existing pinned native theorem citations; F00/E and F00GHI finite audits rerun |

The declared construction choices are visible in [CORE_LEDGER.json](certificate/CORE_LEDGER.json): the typed recursive source, logical/set-theoretic proof language, signed role completion, coefficient-matching pairing, full finite-prefix observation, and halving as metric normalization. Theorems prove consequences of these definitions. They do not prove these choices uniquely follow from the original prose, nor select a physical metric or clock. Countable scalar and history completions are constructed explicitly; no hidden physical probability or noise law is admitted.

Exact checks cover 8 derived signed operators, all 512 group triples, all 32 pairing basis cases, 16 cut-corner basis cases, 13 symbolic proof replays, 63 finite histories and all 250,047 ultrametric triples, 270 potential assignments on four graph types, and the derived 8-vertex/16-edge role graph with 9 independent comparison cycles. Incorrect closure, history erasure, and altered proof certificates are rejected. General infinite/analytic statements have written proofs or named source theorems; finite tests do not prove them. No Lean/Coq proof or external peer review is claimed.

Reproduce against the pinned canonical Recognition Kernel checkout at `3cc5a33b05c16d59c90994ddda69dedc0d392424`:

```bash
python3.12 certificate/verify.py --rkf-root /absolute/path/to/Recognition-Kernel-Framework
```

The verifier checks local and canonical source hashes before loading the adapters, recomputes the witnesses, replays the unchanged RKF engine, and reruns both native analytic audit packets without changing upstream files. It writes `certificate/VERIFICATION.json` and `certificate/NATIVE_REPLAY.json`. Python exact integers/rationals encode proved cut counts; engine complex-like multiplication is used only after its role-derived law is established. The role derivation itself uses only radial coefficients, which the bridge checks explicitly.

To rebuild the readable manuscript, run `pdflatex -interaction=nonstopmode -halt-on-error emk_topology_foundation.tex` twice. The shipped PDF is bound to the packet by its hash; a new TeX build can have different PDF metadata. Rebuilding the PDF is separate from verifying the shipped packet.

Sources are cited by exact commit, path, blob and SHA-256 in [SOURCE_PINS.json](certificate/SOURCE_PINS.json). MP NR-7/NR-8 receive credit for cut refinements/completion, and NR-9 for the scalar-shadow quotient under its declared UGD law. Their private source bodies are not republished. NR-7 originally defines an oriented action; that fact is stated in the citation, and C4 supplies the role construction here. RKF F00/E, F00G and RH F00J/F00K are native analytic sources with their completion premises stated. MP NR-13's later Haar, continuous `L¹/L²` and Young-type machinery is labeled **analytic import, excluded from this core** in its citation.

The original source bytes are preserved in [source/emk_topology.original.tex](source/emk_topology.original.tex). [ORIGINAL_CLAIM_LEDGER.json](certificate/ORIGINAL_CLAIM_LEDGER.json) keeps every original title, location, content hash, audit reason and corrected disposition. A withdrawn physical or infinite-series claim remains withdrawn; it is never counted as proved merely because the new core is closed.
