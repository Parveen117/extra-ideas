# Way-4: Operator evolution

Source: the section of that name in [4ways.tex](../4ways.tex), with later shared sections on canonical operators and unified evolution. [source-excerpt.tex](source-excerpt.tex) preserves the principal section verbatim as a LaTeX fragment.

## Core idea

Describe states through the transformations that act on them: composition, generators, flows, commutators, and spectral data.

## Research judgment

This is the strongest immediate working machinery in the draft. Operators can encode order-sensitive actions, while a bare normalization does not determine such actions. This judgment does not validate the manuscript's current Fredholm, regularization, or completeness claims.

## A precise starting example

On `ell^2(N; C)`, define `S e_i=i e_i` and `L e_i=e_{i+1}`. On an appropriate common domain,

`[S,L]=L`, and consequently `LS=(S-I)L`.

The second identity corrects the sign in the draft. On finitely supported vectors, these relations can be checked directly. The scale operator is unbounded, so domains matter when extending the statements.

## What needs development

Start with an abstract algebra and its actions if the intended foundation is independent of Hilbert space. A complex Hilbert representation then needs its own construction and hypotheses; it cannot be obtained merely by writing `ell^2(K)` for a general valued field.

For evolution, specify whether the generator is bounded, an unbounded semigroup generator, or a formal derivation. For Fredholm equations, distinguish invertibility from singularity and give the compatibility condition for the inhomogeneous problem. For the unified composition, match every map's domain and codomain.

## Next target

Construct one native operator family and determine which spectral or return quantities survive admissible reference changes. Use those quantities to investigate Way-3's selection problem, then derive a Way-2 response. The [lambda note](../03-lambda-reference/IDENTIFIABILITY.md) shows why an adjustable ladder spacing by itself does not determine a fundamental constant.

Role in the four-route program: the main calculation and proof tool. See [the shared assessment](../ASSESSMENT.md) for the other operator corrections.

## R1 development and reproduction

[response_transport.py](response_transport.py) implements the finite scalar return model with exact fractions. The labelled edge operators compose to `T(C)=Lambda_C P_i`; their reference changes are diagonal conjugations. No Hilbert-space primitive is used.

Run the focused checks from the repository root with Python 3.11 or 3.12:

```bash
python3.12 -B 04-operator-evolution/verify_r1.py
```

[test_response_transport.py](test_response_transport.py) checks return invariance, reconstruction, operator composition, multiple cycles, conservation rates, and counterexamples. [R1_VERIFICATION.json](R1_VERIFICATION.json) records the result. The general proofs and scope are in [RETURN_INVARIANTS_R1.md](../03-lambda-reference/RETURN_INVARIANTS_R1.md).
