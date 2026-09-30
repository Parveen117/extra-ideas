# Extra Ideas

Four research routes extracted from `4ways.tex`, organized for further development.

Research owner: **Monty Dabas**. Initial organization and assessment: 30 September 2026.

Public repository: [Parveen117/extra-ideas](https://github.com/Parveen117/extra-ideas).

## The four folders

| Folder | Route | Main question |
| --- | --- | --- |
| [01-intrinsic-numbers](01-intrinsic-numbers/README.md) | Way-1: intrinsic multidimensional numbers | What algebra carries the states and their components? |
| [02-relational-response](02-relational-response/README.md) | Way-2: relational and response numbers | Which changes, constraints, and seam relations are invariant? |
| [03-lambda-reference](03-lambda-reference/README.md) | Way-3: a constant reference lambda | Can the reference be selected by the structure, and yield an identifiable prediction? |
| [04-operator-evolution](04-operator-evolution/README.md) | Way-4: operator evolution | What transformations, compositions, and spectra does the structure support? |

Each folder contains a route guide and a verbatim LaTeX excerpt of its principal section. The excerpts are fragments, not standalone papers. Shared foundations and applications remain in the complete [original manuscript](4ways.tex).

## Latest development: R4 native source bridge and closure checks

[R4: The native source of R3's bond](03-lambda-reference/NATIVE_BOND_BRIDGE_R4.md) connects R3 to the existing RKF paired-depth response. At the source's exact cut, `F=I+KR` and `B=F/2`. With the existing grading `A=K`, the removed-space rule follows: `ABA=I-B`, and `[A,B]=R` recovers the existing quarter-turn. A coordinate map identifies the full R3 bond and relative metric with this source realization.

This resolves R3's separate complement choice **within the supplied exact-cut source family**. It does not select that family, its aperture, or its remaining coupling universally. The quarter-turn and native cut synthesis are earlier results, not new R4 discoveries.

R4 adds a linked raw-defect check, `J_x^2+I=-4(B_x^2-B_x)`, and transfers the native solver's finite-aperture enclosures to both defects. Actual source examples show why normalizing a commutator or observing one exactly closed finite aperture can give a misleading closure verdict.

The 12 focused integration tests pass, including 16 native polynomial proof replays. Reproduce with a separate checkout of RKF commit `3cc5a33b05c16d59c90994ddda69dedc0d392424`:

```bash
python3.12 -B 04-operator-evolution/verify_r4.py --rkf-root ../rkf-r4
```

The verifier checks the four upstream runtime hashes before executing the unchanged solver. See the [source pins](04-operator-evolution/R4_SOURCE_PINS.json), [verification record](04-operator-evolution/R4_VERIFICATION.json), and [cross-repository lineage assessment](CROSS_REPO_LINEAGE.md). The assessment identifies R2's compression identity as an existing theorem's special case and distinguishes R3's representation from earlier iota derivations.

## R3: seam-selected bond and complex structure

[R3: A typed seam, its bond, and a linear complex structure](03-lambda-reference/SEAM_BOND_COMPLEX_STRUCTURE_R3.md) supplies a conditional selection rule for R2's bond: retain the seam `L` and remove its Aghora image `A L`, when they are complementary. This **added rule** uniquely determines `B` and proves:

- `ABA=I-B` and `K_Gamma=[A,B]` satisfies `K_Gamma^2=-I`.
- In a real two-mode realization, requiring `A` to be an isometry and `B` orthogonal fixes a positive relative metric up to overall scale.
- Requiring the Aghora-odd generator to preserve that metric gives `G=omega K_Gamma`, with `omega` still free.
- The exponential return compresses to `BR_tB=sin(omega t)B`; no-leakage returns allow either sign.

These conditions construct a standard linear complex structure from the specified seam and involution. They are not implied by the source ratio alone. The metric is on a two-mode state space; no physical metric, absolute rate, or universal `lambda_*` has been selected.

Read the [complex-coordinate action](01-intrinsic-numbers/SEAM_COMPLEX_COORDINATES_R3.md), [typed response](02-relational-response/TYPED_SEAM_R3.md), and [exact implementation](04-operator-evolution/seam_bond.py). Reproduce the 18 focused checks with:

```bash
python3.12 -B 04-operator-evolution/verify_r3.py
```

[R3_VERIFICATION.json](04-operator-evolution/R3_VERIFICATION.json) records exact examples and counterexamples, hashes the implementation and proof, and checks the preserved R1/R2 evidence. Its chosen rational example has compressed coefficient `3/5` and removed squared-norm fraction `16/25`. A counterexample demonstrates that zero scalar defect can hide leakage if the generator lacks metric compatibility.

## R2: Aghora returns and projection

[R2: Aghora return law and the effect of a projection](03-lambda-reference/AGHORA_RETURN_R2.md) uses the manuscript's `A^2=I`, `AGA=-G`, and bond projection. With the stated specialization `U_t=exp(tG)` and constructed protocol `R_t=A U_t`, it proves:

- `R_t^2=I`; the return is similar to `A`, so its eigenvalues remain in `{+1,-1}`.
- A nonzero odd generator requires both Aghora sectors and at least two modes.
- A bond selects one sign only when its retained space lies in the corresponding return eigenspace.
- The projected return satisfies the exact identity `B-(BR_tB)^2=BR_t(I-B)R_tB`.

Under the additional orthogonality and adjoint assumptions stated in the proof, the last expression is a nonnegative squared-norm defect. A rational example has full-return signs `+1,-1`, a projected coefficient `4/5`, and omitted squared norm `9/25`. Those fractions are illustrative, not universal constants.

Read the [mode construction](01-intrinsic-numbers/AGHORA_MODES_R2.md), [projected response](02-relational-response/PROJECTED_RESPONSE_R2.md), and [exact matrix implementation](04-operator-evolution/aghora_return.py). Reproduce its 15 focused checks with:

```bash
python3.12 -B 04-operator-evolution/verify_r2.py
```

[R2_VERIFICATION.json](04-operator-evolution/R2_VERIFICATION.json) records the checked source hashes. R2 supplies a particular constrained return family. R3 adds an explicit conditional bond rule, and R4 realizes that rule through an existing native exact-cut response. Universal source selection and a value of `lambda_*` remain unresolved. The constructed return is not yet identified with the source's full six-stage closure.

## R1: reference-independent return invariants

[R1: Reference-independent return multipliers](03-lambda-reference/RETURN_INVARIANTS_R1.md) constructs a finite scalar transport model joining the four routes. It proves:

- Closed-return products are invariant under independent local reference rescalings.
- A spanning tree leaves exactly `m-n+1` nonzero scalar return parameters, which completely classify that model up to the stated reference changes.
- Such an invariant is constant during evolution precisely when its cycle logarithmic rate vanishes; all returns are constant when edge rates are vertex differences.
- Labelled coordinate operators reproduce every response path and its return.

The [coordinate model](01-intrinsic-numbers/COORDINATE_TRANSPORT_R1.md), [response interpretation](02-relational-response/RESPONSE_RETURNS_R1.md), and [exact implementation](04-operator-evolution/response_transport.py) make the correspondence explicit. Run the focused checks with Python 3.11 or 3.12:

```bash
python3.12 -B 04-operator-evolution/verify_r1.py
```

See [R1_VERIFICATION.json](04-operator-evolution/R1_VERIFICATION.json) for the recorded result and source hashes. The model assumes reciprocal scalar transport and independently changeable local references. Its return values are computable from edge data; a universal value of `lambda_*` remains unselected. R2 and R3 investigate additional operator and seam structure. The proof document acknowledges the established gain-graph mathematics underlying this construction.

## Research judgment

**Way-4 currently offers the strongest working machinery. Way-3 offers an ambitious foundational question, but its constant has not yet been determined by the draft.** Way-2 is the bridge to response predictions; Way-1 supplies the coordinate algebra.

The proposed development priority is to investigate Way-3 using Way-4, then express the resulting invariant in Way-2. In particular, distinguish a fixed reference from a uniquely selected constant. With unrestricted coefficients, writing `X_i = n_i lambda_*^{k_i}` alone does not select `lambda_*`.

See [ASSESSMENT.md](ASSESSMENT.md) for the reasoning and corrections, and [the lambda identifiability note](03-lambda-reference/IDENTIFIABILITY.md) for an exact rescaling argument and a concrete next research target.

## Status of claims

This is a research workspace. The source manuscript is preserved byte for byte, including its original theorem labels, wording, and placeholder bibliography. Preservation is not an endorsement of every claim. Several statements require additional hypotheses, and some have explicit counterexamples.

The unrestricted four-way equivalence remains a proposal. R1 now supplies an explicit correspondence for a specified finite scalar transport model, with normalization retaining its reference class. No value of a physical constant or experimentally validated new physical law is established.

## Selection milestone

R1 specifies a scalar model, its reference changes, and its complete return invariants. R2 instantiates the existing cut-graded return and compression calculus. R3 supplies an explicit conditional seam/bond realization. R4 connects it to an existing native exact-cut response and checks closure without discarding its amplitude or aperture dependence. Next supply a source-selection or dynamical law that distinguishes the remaining profiles, apertures, rates or phases. Invariance, conservation, and universal value selection remain separate claims.

The name `lambda` currently serves several roles. Use `lambda_*` for the fixed reference, `chi_XY` for a response derivative, `gamma` for a decay rate, and `epsilon` for observational resolution until a theorem relates them.

## Provenance

[MANIFEST.json](MANIFEST.json) records the original source checksum and the line ranges of the excerpts. R4 executes a pinned upstream implementation from a separate checkout; its source files are referenced rather than duplicated here. [CROSS_REPO_LINEAGE.md](CROSS_REPO_LINEAGE.md) records inherited results and the limited contribution of each development.
