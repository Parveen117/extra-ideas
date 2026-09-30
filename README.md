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

## Latest development: R1 return invariants

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

R1 specifies a model, its reference changes, and its complete return invariants. Next derive a native return law that fixes a distinguished invariant across admissible systems. Give its existence and uniqueness proof and a response prediction that could distinguish it from competing choices. Invariance under reference changes alone does not imply a universal or conserved value.

The name `lambda` currently serves several roles. Use `lambda_*` for the fixed reference, `chi_XY` for a response derivative, `gamma` for a decay rate, and `epsilon` for observational resolution until a theorem relates them.

## Provenance

[MANIFEST.json](MANIFEST.json) records the original source checksum and the line ranges of the excerpts. No external repository material has been imported.
