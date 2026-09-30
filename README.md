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

## Latest development: R2 Aghora returns and projection

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

[R2_VERIFICATION.json](04-operator-evolution/R2_VERIFICATION.json) records the checked source hashes. R2 supplies a particular constrained return family; the native bond-selection law and a universal value of `lambda_*` remain unresolved. The constructed return is not yet identified with the source's full six-stage closure.

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

See [R1_VERIFICATION.json](04-operator-evolution/R1_VERIFICATION.json) for the recorded result and source hashes. The model assumes reciprocal scalar transport and independently changeable local references. Its return values are computable from edge data; a universal value of `lambda_*` remains unselected. The native return operation and a law selecting its value are the next research target. The proof document acknowledges the established gain-graph mathematics underlying this construction.

## Research judgment

**Way-4 currently offers the strongest working machinery. Way-3 offers an ambitious foundational question, but its constant has not yet been determined by the draft.** Way-2 is the bridge to response predictions; Way-1 supplies the coordinate algebra.

The proposed development priority is to investigate Way-3 using Way-4, then express the resulting invariant in Way-2. In particular, distinguish a fixed reference from a uniquely selected constant. With unrestricted coefficients, writing `X_i = n_i lambda_*^{k_i}` alone does not select `lambda_*`.

See [ASSESSMENT.md](ASSESSMENT.md) for the reasoning and corrections, and [the lambda identifiability note](03-lambda-reference/IDENTIFIABILITY.md) for an exact rescaling argument and a concrete next research target.

## Status of claims

This is a research workspace. The source manuscript is preserved byte for byte, including its original theorem labels, wording, and placeholder bibliography. Preservation is not an endorsement of every claim. Several statements require additional hypotheses, and some have explicit counterexamples.

The unrestricted four-way equivalence remains a proposal. R1 now supplies an explicit correspondence for a specified finite scalar transport model, with normalization retaining its reference class. No value of a physical constant or experimentally validated new physical law is established.

## Selection milestone

R1 specifies a scalar model, its reference changes, and its complete return invariants. R2 constructs an involutive operator-return family from the source relations and identifies how a bond can change its observed response. Next derive the bond or seam selection law: determine whether it selects a sign sector or fixes a specific noncommuting projection. Invariance, conservation, and universal value selection remain separate claims.

The name `lambda` currently serves several roles. Use `lambda_*` for the fixed reference, `chi_XY` for a response derivative, `gamma` for a decay rate, and `epsilon` for observational resolution until a theorem relates them.

## Provenance

[MANIFEST.json](MANIFEST.json) records the original source checksum and the line ranges of the excerpts. No external repository material has been imported.
