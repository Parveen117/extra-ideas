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

Extend the [R7 tensor foundation](../02-relational-response/EMK_TENSOR_CALCULUS_R7.md) with a native direction/soldering rule, a metric policy and higher exterior operations before physical applications. KIR transport is already lifted to explicit tensor types; its finite group action and infinitesimal Lie action must remain distinct. The [lambda note](../03-lambda-reference/IDENTIFIABILITY.md) retains the separate selection obligation.

Role in the four-route program: the main calculation and proof tool. See [the shared assessment](../ASSESSMENT.md) for the other operator corrections.

## R1 development and reproduction

[response_transport.py](response_transport.py) implements the finite scalar return model with exact fractions. The labelled edge operators compose to `T(C)=Lambda_C P_i`; their reference changes are diagonal conjugations. No Hilbert-space primitive is used.

Run the focused checks from the repository root with Python 3.11 or 3.12:

```bash
python3.12 -B 04-operator-evolution/verify_r1.py
```

[test_response_transport.py](test_response_transport.py) checks return invariance, reconstruction, operator composition, multiple cycles, conservation rates, and counterexamples. [R1_VERIFICATION.json](R1_VERIFICATION.json) records the result. The general proofs and scope are in [RETURN_INVARIANTS_R1.md](../03-lambda-reference/RETURN_INVARIANTS_R1.md).

## R2 development and reproduction

[aghora_return.py](aghora_return.py) supplies exact rational matrices for the Aghora return and bond projection. It evaluates nilpotent exponentials and a separately identified rational reversible family. No approximate matrix exponential is used.

```bash
python3.12 -B 04-operator-evolution/verify_r2.py
```

The 15 [focused checks](test_aghora_return.py) cover the return identities, both sign sectors, projection defects, reference changes, and counterexamples when hypotheses are missing. [R2_VERIFICATION.json](R2_VERIFICATION.json) records the run. Read [the proof](../03-lambda-reference/AGHORA_RETURN_R2.md) for the distinction between algebraic identities, optional norm statements, and the unresolved native selection law.

## R3 development and reproduction

[seam_bond.py](seam_bond.py) constructs the projection onto a two-mode seam along its Aghora image, the resulting complex structure, a compatible relative metric, and exact rational returns. The rule for the removed subspace and the metric conditions are explicit additions to the source.

```bash
python3.12 -B 04-operator-evolution/verify_r3.py
```

The 18 [focused checks](test_seam_bond.py) cover the derived operators, changes of coordinates, complex multiplication, the response and norm balance, both return signs, and counterexamples without the added hypotheses. [R3_VERIFICATION.json](R3_VERIFICATION.json) hashes R3 and its imported R2 implementation, and checks the preserved earlier verification hashes. The [proof](../03-lambda-reference/SEAM_BOND_COMPLEX_STRUCTURE_R3.md) distinguishes the algebraic construction from optional metric and dynamical requirements. No floating-point approximation is used.

## R4 source bridge and reproduction

[native_bond_bridge.py](native_bond_bridge.py) connects the raw response to bond/commutator defects and finite error budgets. [native_bond_probe.cjs](native_bond_probe.cjs) calls the unchanged native solver and proof replayer in a separate RKF checkout. The verifier refuses altered upstream bytes before execution.

Use Python 3.11 or 3.12 and Node.js, with RKF at the commit pinned in [R4_SOURCE_PINS.json](R4_SOURCE_PINS.json):

```bash
python3.12 -B 04-operator-evolution/verify_r4.py --rkf-root ../rkf-r4
```

The 12 [integration tests](test_native_bond_bridge.py) include 16 native polynomial replays, the actual Schur witness, the R3 coordinate/metric map, finite enclosures, and false-positive closure controls. [R4_VERIFICATION.json](R4_VERIFICATION.json) records the source hashes, results and preserved R1/R2/R3 evidence. Read the [proof](../03-lambda-reference/NATIVE_BOND_BRIDGE_R4.md) and [lineage assessment](../CROSS_REPO_LINEAGE.md) for what is inherited and what the bridge adds.

## R5 profile recovery and reproduction

[profile_recovery.py](profile_recovery.py) recovers a positive depth prefix by reciprocal-series stripping and reports that its tail remains unknown. Inputs are exact coefficients of `z, z^3, ...` for a calibrated common multiplier of paired products. A finite rational-polynomial companion is checked against the original source, while [native_profile_probe.cjs](native_profile_probe.cjs) reuses its forward responses, synthesis and interval routines unchanged.

```bash
python3.12 -B 04-operator-evolution/profile_recovery.py '[3, -18, 216]'
python3.12 -B 04-operator-evolution/verify_r5.py --rkf-root ../rkf-r4
```

The CLI returns `[3,2,3]` and explicitly leaves infinite periodicity unproved. The 14 [focused tests](test_profile_recovery.py) include heterogeneous reconstruction, native forward checks, the equal-cut ambiguity witness, disjoint held-out response intervals, gain ambiguity, and finite tail budgets. [R5_VERIFICATION.json](R5_VERIFICATION.json) records the results and preserved R1–R4 evidence. See the [proof](../03-lambda-reference/PROFILE_RECOVERY_R5.md) and [source pins](R5_SOURCE_PINS.json) for assumptions and classical lineage. No upstream runtime is copied into this repository.

## R6 interval recovery and reproduction

[interval_profile_recovery.py](interval_profile_recovery.py) propagates supplied rational coefficient intervals through R5's inverse, returns certified coupling ranges and a stopping status, and predicts a response with both prefix error and the remaining tail. R5's native adapter and all earlier proof/code files stay byte-identical.

```bash
python3.12 -B 04-operator-evolution/verify_r6.py --rkf-root ../rkf-r4
```

The 10 [focused checks](test_interval_profile_recovery.py) include eight prediction extrema replayed through the unchanged native solver, error-aware source discrimination and refusal cases. See the [proof](../03-lambda-reference/INTERVAL_PROFILE_RECOVERY_R6.md) and [verification record](R6_VERIFICATION.json). Use Python 3.11 or 3.12 and the pinned R4 upstream runtime.

## R7 tensor foundation and reproduction

[emk_tensor_calculus.py](emk_tensor_calculus.py) implements exact typed tensors, slot transport, generator derivatives, contraction, symmetry/wedge, metric conversion, seam blocks through projectors, finite increments, ordered integration, curvature, invariant symmetric forms and carrier/sheet return audits. Its scalar field is Q; the KIR operators and integer winding register are separately typed.

```bash
python3.12 -B 04-operator-evolution/verify_r7.py --publications-root ../Publications
```

Use Python 3.11 or 3.12 with a separate Publications checkout at the commit in [R7_SOURCE_PINS.json](R7_SOURCE_PINS.json). [verify_r7.py](verify_r7.py) rejects modified EMK2/EMKT1/EMKT2 source bytes before importing them. The 16 [focused tests](test_emk_tensor_calculus.py) include 30 rank/variance signatures, metric blindness, the fixed-metric obstruction, curved D², finite Bianchi and independent sheet closure. [R7_VERIFICATION.json](R7_VERIFICATION.json) records the result and all earlier recorded hashes. The [proof](../02-relational-response/EMK_TENSOR_CALCULUS_R7.md) separates the executable finite layer from optional smooth and physical structures.
