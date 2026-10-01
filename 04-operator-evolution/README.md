# Way-4: Operator evolution

Source: the section of that name in [4ways.tex](../4ways.tex), with later shared sections on canonical operators and unified evolution. [source-excerpt.tex](source-excerpt.tex) preserves the principal section verbatim as a LaTeX fragment.

## R15 cut/operator and topology verification

[R15's proof and source audit](../02-relational-response/EMK_TOPOLOGY_R15.md) accompany [emk_topology_core.py](emk_topology_core.py), the [unchanged-engine caller](r15_native_topology_probe.cjs), [48-claim ledger](certificates/r15/CLAIM_LEDGER.json) and [source pins](certificates/r15/SOURCE_PINS.json). The full original manuscript is **not certified as written**; corrected statements have separate written proofs and scoped exact checks.

```bash
python3.12 -B 04-operator-evolution/verify_r15.py \
  --rkf-root /path/to/Recognition-Kernel-Framework \
  --publications-root /path/to/Publications
```

[R15_VERIFICATION.json](R15_VERIFICATION.json) records exhaustive checks on 3,678 finite systems/observations, 42 signed phase basis vectors, ten native symbolic replays, negative controls and an exact unchanged UGD-1 certificate replay. No upstream runtime is copied and no old upstream pin is overwritten.

## R14 executable scalar specialization

[R14's foundational interpretation is superseded](../02-relational-response/R14_SCOPE_CORRECTION.md). Its exact scalar calculations remain reproducible; they do not constitute the full generalized UGD source.

[R14](../02-relational-response/NATIVE_SOURCE_FOUNDATION_R14.md) and its [citation ledger](../02-relational-response/DERIVATION_SOURCE_LEDGER_R14.md) link the existing RKF/RH foundation to exact source-driven memory, record-event content, covariance and phase geometry. The adapter [native_source_foundation.cjs](native_source_foundation.cjs) calls the canonical RKF arithmetic source unchanged.

```bash
node 04-operator-evolution/verify_r14.cjs \
  --rkf-root /path/to/Recognition-Kernel-Framework \
  --output /tmp/R14_VERIFICATION.json
```

The [certificate](R14_VERIFICATION.json) records 30 checks, 36 histories, 28 hidden recoveries, 54 complete record trees, 36 metric tangents, nine observer completions and eight rejected mathematical mutations. [Pins](R14_SOURCE_PINS.json) identify the sources. [Upstream replay status](R14_UPSTREAM_REPLAYS.json) includes an old combined-engine certificate hash failure; no full upstream-engine PASS or unconditional physical foundation is claimed.

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

Follow [R13](../02-relational-response/COMPLEMENT_MEMORY_NOISE_R13.md): complementary transport now generates memory and endogenous covariance, and native norm-weighted binary records realize the observer/metric equation. Next match the preparation, monitoring policy and record coupling to an actual native source and test its colored-force, event and recovery predictions. Preserve signed information and the complete R10 connection contract.

## R13 complementary memory/noise certificates

[complement_memory_noise.py](complement_memory_noise.py) implements exact block elimination, complementary covariance, balanced unit preparations, hidden-state recovery, native first-exit weights, event Fisher information and a non-Gaussian binary realization of the completed seam metric. [r13_native_complement_probe.cjs](r13_native_complement_probe.cjs) calls the unchanged KIR presentation and N03 engine for 11 algebra replays and three observer completions.

~~~bash
python3.12 -B 04-operator-evolution/verify_r13.py \
  --publications-root ../Publications \
  --rkf-root ../Recognition-Kernel-Framework
~~~

[R13_VERIFICATION.json](R13_VERIFICATION.json) binds 31 exact tests, six rejected mathematical mutations, three rejected native alterations and unchanged R1–R12/master evidence. See the [proof](../02-relational-response/COMPLEMENT_MEMORY_NOISE_R13.md), [pins](R13_SOURCE_PINS.json) and [native packet](R13_NATIVE_CERTIFICATE.json). Native binary likelihood replaces the need for Gaussian noise in this construction; coherent retention and monitored reset remain distinct protocols.

## R12 observer/metric selection certificates

[observer_metric_foundation.py](observer_metric_foundation.py) derives the tagged Fisher form, future observer quotient, descended transports, tag-erasure ledger, transport-derived metric two-jet and quadratic seam geometry. [r12_native_observer_metric_probe.cjs](r12_native_observer_metric_probe.cjs) calls unchanged RKF routines for nine algebra replays and four N03 observer completions.

~~~bash
python3.12 -B 04-operator-evolution/verify_r12.py \
  --publications-root ../Publications \
  --rkf-root ../Recognition-Kernel-Framework
~~~

[R12_VERIFICATION.json](R12_VERIFICATION.json) records 28 tests, six rejected mathematical mutations, three rejected native alterations and preserved R1–R11/master hashes. See the [proof](../02-relational-response/OBSERVER_METRIC_FOUNDATION_R12.md), [pins](R12_SOURCE_PINS.json) and [native packet](R12_NATIVE_CERTIFICATE.json). The operational metric is derived within its supplied event/noise law; universal physical metric and clock selection remain open.

## R11 spectral curvature observer certificates

[spectral_curvature_observer.py](spectral_curvature_observer.py) implements the admitted curvature family, exact finite loops, characteristic and marked determinants, minimum target repair, catalogue search, cubic connected coefficients, conditional response geometry and noise decisions. [r11_native_spectral_probe.cjs](r11_native_spectral_probe.cjs) calls the unchanged RKF engine for eight equalities and one distinctness replay.

    python3.12 -B 04-operator-evolution/verify_r11.py \
      --publications-root ../Publications \
      --rkf-root ../Recognition-Kernel-Framework

Use Python 3.11 or 3.12, Node and [R11_SOURCE_PINS.json](R11_SOURCE_PINS.json). The verifier binds the original spectral PDFs and checks native bytes before replay. [R11_VERIFICATION.json](R11_VERIFICATION.json) records 31 tests, nine replays, six rejected math mutations, two rejected native alterations and preserved R1–R10/master evidence. See the [proof](../02-relational-response/SPECTRAL_CURVATURE_OBSERVER_R11.md) and [native packet](R11_NATIVE_CERTIFICATE.json). Completeness is model-relative; physical implementation and erased history are separate questions.

## R10 native curvature descent certificates

[native_curvature_descent.py](native_curvature_descent.py) implements constant and smooth-jet sourced Bianchi, metric two-jet Levi-Civita construction, surjective typed observer descent, Riemann/distortion/excursion decomposition and finite retained-memory witnesses. [r10_native_curvature_probe.cjs](r10_native_curvature_probe.cjs) calls the unchanged canonical engine to replay three single-carrier and two typed quotient proofs; no engine is copied here.

```bash
python3.12 -B 04-operator-evolution/verify_r10.py \
  --publications-root ../Publications \
  --rkf-root ../Recognition-Kernel-Framework
```

Use Python 3.11 or 3.12, Node and the commits in [R10_SOURCE_PINS.json](R10_SOURCE_PINS.json). [R10_VERIFICATION.json](R10_VERIFICATION.json) records 25 exact tests, five symbolic replays, six rejected mathematical mutations, two rejected native alterations and unchanged R1–R9/master evidence. See the [proof](../02-relational-response/NATIVE_CURVATURE_DESCENT_R10.md) and [native packet](R10_NATIVE_CERTIFICATE.json) for the local-jet, constant-cut and conditional tangent-sector scope.

## R9 native balance and information certificates

[emk_curvature_balance.py](emk_curvature_balance.py) implements cut grading, recognition events, curvature-sector selection, a mean/variance budget, whole-family information pushforward and finite measurement ledgers. [r9_native_balance_probe.cjs](r9_native_balance_probe.cjs) calls the unchanged **RKF/operator_foundation** engine for six general involution identities and three KIR certificates. The derived cut and information operations are consumed from unchanged EMK-T1/QTH-1 source.

```bash
python3.12 -B 04-operator-evolution/verify_r9.py \
  --publications-root ../Publications \
  --rkf-root ../Recognition-Kernel-Framework
```

Use Python 3.11 or 3.12, Node, and the separate checkouts at commits in [R9_SOURCE_PINS.json](R9_SOURCE_PINS.json). [R9_VERIFICATION.json](R9_VERIFICATION.json) records 24 passing exact tests, nine symbolic replays, six rejected math mutations, two rejected native alterations and preserved R1–R8/master-review evidence. [R9_NATIVE_CERTIFICATE.json](R9_NATIVE_CERTIFICATE.json) contains the pinned symbolic contracts and replayable witnesses. The [proof](../02-relational-response/CURVATURE_BALANCE_R9.md) states the exchange-symmetry hypothesis and keeps recognized curvature, information cost, branch recovery and Riemann curvature distinct.

## R8 exact curvature and observation adapter

[emk_curvature_observation.py](emk_curvature_observation.py) reuses the R2 matrix arithmetic and R7 chart for compression, odd-target reconstruction and connection jets. [r8_native_curvature_probe.cjs](r8_native_curvature_probe.cjs) calls the unchanged canonical **RKF/operator_foundation** engine; no engine implementation is copied here. It replays a symbolic general compression proof, two KIR identities and a rank-2 to rank-4 future observer completion.

```bash
python3.12 -B 04-operator-evolution/verify_r8.py --rkf-root ../Recognition-Kernel-Framework
```

Use the separate RKF checkout at the commit in [R8_SOURCE_PINS.json](R8_SOURCE_PINS.json), Python 3.11 or 3.12, and Node. [R8_VERIFICATION.json](R8_VERIFICATION.json) records 23 passing exact tests, six rejected mathematical mutations, three rejected altered native contracts/results, and preservation of earlier evidence. [R8_NATIVE_CERTIFICATE.json](R8_NATIVE_CERTIFICATE.json) contains the replayable native inputs, rewrite witnesses and observer result. General written proofs and the exact computation have explicit assumptions; physical validation and new Lean formalization remain separate evidence levels.

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
