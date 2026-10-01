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

## R14: native source foundation and cross-repository derivation ledger

[R14](02-relational-response/NATIVE_SOURCE_FOUNDATION_R14.md) derives the radial/complement split from native dagger and the transport from multiplication by `z=a+iota b`. Existing RH native period/polar/logarithm results derive theta from that source. The same source gives exact memory, consistent record-event content, and a phase metric; source-cell counts give unresolved-force covariance without imposing a balanced or Gaussian preparation.

The native normalized-amplitude metric and the separately derived record-logarithm Hessian agree up to the proved factor four:

\[
h=\frac{(a\,db-b\,da)^2}{(a^2+b^2)^2}=d\theta^2,
\qquad g_{\mathrm{record}}=4h.
\]

The [derivation ledger](02-relational-response/DERIVATION_SOURCE_LEDGER_R14.md) integrates **24 commit-pinned references** from RKF and RH, including mathematics developed for Yang–Mills. Each citation states its native inputs or classical imports and its actual evidence level. The current source-selection and physical event-actualization questions remain explicit; R14 does not relabel a source tuple as an unconditional derivation.

[Thirty exact checks and eight mathematical mutation controls pass](04-operator-evolution/R14_VERIFICATION.json). The [upstream replay record](04-operator-evolution/R14_UPSTREAM_REPLAYS.json) includes passing F00/E, GHI and 45-test RH foundation replays, and an explicitly failed old combined-engine certificate hash check. R1–R13 evidence is unchanged.

```bash
node 04-operator-evolution/verify_r14.cjs \
  --rkf-root /path/to/Recognition-Kernel-Framework \
  --output /tmp/R14_VERIFICATION.json
```

## R13: complementary memory derives noise and event information

[R13](02-relational-response/COMPLEMENT_MEMORY_NOISE_R13.md) constructs the statistical ingredients that R12 supplied. Exact hidden-state elimination generates memory and an unresolved complementary force. A balanced unit-norm native preparation derives its colored covariance. Additive native norm readout derives continuation and exit weights \(r=\cos^2\theta,\ q=\sin^2\theta\), with a monitored first-exit law.

Native balanced binary meters then realize the completed observer/metric selection equation **without Gaussian noise or external record precision**. Matching the sheet gain to the complementary sine leg gives

\[
g=(1+\cos^2\theta\,v^2)\,du^2+dv^2.
\]

The cut/transport commutator, complementary norm loss and returning memory all use the same sine amplitude. Coherent retention, monitored reset, hidden-force fluctuations and binary event fluctuations have explicit separate contracts.

The [verification record](04-operator-evolution/R13_VERIFICATION.json) reports **31 exact tests**, 11 native replays, three canonical N03 completions, six rejected mathematical mutations and three rejected native alterations. All R1–R12/master evidence is preserved. See the [pins](04-operator-evolution/R13_SOURCE_PINS.json) and [native packet](04-operator-evolution/R13_NATIVE_CERTIFICATE.json).

~~~bash
python3.12 -B 04-operator-evolution/verify_r13.py \
  --publications-root ../Publications \
  --rkf-root ../Recognition-Kernel-Framework
~~~

The derivation selects noise, event weights and information geometry within a concrete native construction. Preparation, rotational phase, norm-readout policy and matched channel identification are specified; their physical selection is a further question.


## R12: deriving the observer and metric

[R12](02-relational-response/OBSERVER_METRIC_FOUNDATION_R12.md) derives a minimum future observer and a unique operational metric from one native event experiment. With calibrated terminal noise and a stable stopping law, its tagged Fisher form solves

\[
G=(1-r)D^TWD+r\sum_a p_aT_a^TGT_a,\qquad
\ker G=\bigcap_w\ker(DT_w).
\]

The balanced nilpotent sheet family selects the existing positive quadratic seam geometry:

\[
g=(1+\kappa v^2)\,du^2+dv^2,\qquad \kappa=c^2r/(1-r).
\]

Metric first and second derivatives are derived from the same equation and match the unchanged EMK-G1 geometry. Balanced linear curvature can vanish while the second-moment metric remains curved. Actual event-tag erasure loses that warp at the common zero state, with an exact positive information ledger; this is a conditional observation law.

The [verification record](04-operator-evolution/R12_VERIFICATION.json) reports **28 exact tests**, nine native algebra replays, four canonical observer completions, six rejected mathematical mutations and three rejected native alterations. R1–R11/master evidence is preserved. See the [pins](04-operator-evolution/R12_SOURCE_PINS.json) and [native certificate](04-operator-evolution/R12_NATIVE_CERTIFICATE.json).

~~~bash
python3.12 -B 04-operator-evolution/verify_r12.py \
  --publications-root ../Publications \
  --rkf-root ../Recognition-Kernel-Framework
~~~

This is a foundation derivation under a specified likelihood, event law and tangent chart. Their physical selection remains open; information distance does not automatically select spacetime signature or force the original native connection to be Levi–Civita.


## R11: spectral blindness and the curvature observer

[R11](02-relational-response/SPECTRAL_CURVATURE_OBSERVER_R11.md) connects both spectral papers in Publications to native curvature and the R10 tangent quotient. It certifies an entire curvature family whose classical geometry is flat and whose characteristic data are constant. Its nonidentity finite order loops can have exactly the spectrum of identity.

Cross-sector insertion markers recover all hidden curvature coordinates. For two visible and h hidden modes, the full declared target needs exactly **2h independent linear scalar responses**. A catalogue without the reverse-sector leg stays blind even with more spectral samples. The spectral paper's cubic connected response and R10's visible excursion use the same closed composition.

The [verification record](04-operator-evolution/R11_VERIFICATION.json) binds **31 passing exact tests**, nine unchanged canonical-engine replays, six rejected mathematical mutations and two rejected native alterations. It includes 81 lower-coupling cases, 405 ordinary and 1,620 marked determinant evaluations, stable noise bounds and model-relative refusal cases. R1–R10/master evidence is preserved. See the [native packet](04-operator-evolution/R11_NATIVE_CERTIFICATE.json) and [source pins](04-operator-evolution/R11_SOURCE_PINS.json).

    python3.12 -B 04-operator-evolution/verify_r11.py \
      --publications-root ../Publications \
      --rkf-root ../Recognition-Kernel-Framework

Use Python 3.11 or 3.12, Node and pinned separate sources. Markers are inserted before compression; erased history requires its own record. Completeness is relative to the declared curvature family, and its response Gram matrix is not automatically a spacetime metric.

## R10: native curvature conservation and the classical tangent sector

[R10](02-relational-response/NATIVE_CURVATURE_DESCENT_R10.md) proves a sourced second Bianchi law for the balance-selected curvature readout: hidden odd couplings account exactly for its visible cyclic current. The full connection retains ordinary Bianchi. A smooth derivative witness and a nonzero three-mode source witness are certified.

Classical Riemann curvature is an exact tangent quotient when a surjective observer intertwines the complete connection jet with an admitted metric's Levi-Civita jet. The visible curvature decomposition is

\[
(F_{ij})_{\mathrm{vis}}=R_{ij}(g)+\mathcal D_{ij}(S)+E_iL_j-E_jL_i.
\]

The certificate rejects pointwise metric/torsion checks and accidental curvature equality as sufficient quotient tests. It also certifies flat Riemann geometry with nonzero hidden native curvature, hidden feedback producing visible curvature, and visible return with retained hidden or sheet memory. The metric jet matches the existing EMK-G1 source at 23 rational positive-domain points; metric selection remains a separate research problem.

The [verification record](04-operator-evolution/R10_VERIFICATION.json) binds **25 passing exact tests**, five unchanged canonical-engine symbolic replays, six rejected mathematical mutations, two rejected native alterations, and preserved R1–R9/master-review evidence. See the [native packet](04-operator-evolution/R10_NATIVE_CERTIFICATE.json) and [source pins](04-operator-evolution/R10_SOURCE_PINS.json).

```bash
python3.12 -B 04-operator-evolution/verify_r10.py \
  --publications-root ../Publications \
  --rkf-root ../Recognition-Kernel-Framework
```

Use Python 3.11 or 3.12, Node and separate pinned source checkouts. These are written proofs, native rewrite certificates and exact finite checks. The canonical engine remains in **RKF/operator_foundation**.

## R9: native curvature balance and the information ledger

[R9](02-relational-response/CURVATURE_BALANCE_R9.md) develops conditional information/curvature balance through the existing derived cut-swap `K`. Within a declared two-branch recognition protocol, invariance under branch exchange uniquely selects equal weights. The resulting balanced observer retains even curvature and removes odd curvature:

\[
\Phi([A,B])=[A_e,B_e]+[A_o,B_o].
\]

For the stated native radial/tangential chart, the `K,R` and `K,RK` curvatures cancel under this observer, while the `R,RK` curvature survives. A clock-free event sequence gives an exact product law for imbalance, including conditional sign reversal and finite balance. The mean-curvature/variance budget records the retained second moment.

The information certificate transforms the complete QTH-1 state family, including its derivatives: `g'=diag(1,lambda^2)` and `B'_xy=lambda^2 m`. Retaining the realized cut branch permits exact inverse recovery; discarding it leaves a quantified positive information loss. Zero mean curvature is also tested against joint parameter-information saturation. Recognition curvature and classical Riemann curvature are documented as different typed objects, with an explicit adapter required to identify them.

The [verification record](04-operator-evolution/R9_VERIFICATION.json) reports **24 passing exact tests**, nine unchanged canonical-engine rewrite replays, six rejected mathematical mutations and two rejected native alterations. R1–R8 and master-review evidence is preserved. See the [symbolic packet](04-operator-evolution/R9_NATIVE_CERTIFICATE.json) and [source pins](04-operator-evolution/R9_SOURCE_PINS.json).

```bash
python3.12 -B 04-operator-evolution/verify_r9.py \
  --publications-root ../Publications \
  --rkf-root ../Recognition-Kernel-Framework
```

Use Python 3.11 or 3.12, Node and the pinned separate checkouts. Exchange symmetry is the endpoint-selection hypothesis; a physical recognition rate is a further model choice. The canonical engine remains in **RKF/operator_foundation**.

## R8: curvature and observation certificates

[R8](02-relational-response/CURVATURE_OBSERVATION_R8.md) links full ordered mismatch, observation, cut sign, dimension restriction and Onsager response flatness through explicit maps. Its central law is

\[
[PAP,PBP]=P[A,B]P-\big(PA(I-P)BP-PB(I-P)AP\big).
\]

The discarded-sector term explains two certified counterexamples: full curvature can be nonzero while reduced generators commute, and commuting full generators can acquire reduced ordered mismatch after projection. On the existing KIR odd-curvature space, two oriented cross-sector scalar channels recover the full current target; the same rows require two more channels for unrestricted future algebra observation. The base-direction pullback rule and constant-response condition `F_ij=L_ji-L_ij` state precisely which flatness claim each observation supports. State-dependent symmetric response and noncommuting flat connections supply the boundary cases.

The unchanged canonical engine remains in **RKF/operator_foundation**. R8 calls it to replay the general idempotent compression identity without a finite carrier, two existing KIR equalities, and a minimum future-observer completion. The [verification record](04-operator-evolution/R8_VERIFICATION.json) reports **23 passing exact tests**, six rejected mathematical mutations, three rejected altered native contracts/results, and preserved R1–R7/master-review hashes. The [complete native packet](04-operator-evolution/R8_NATIVE_CERTIFICATE.json) and [source pins](04-operator-evolution/R8_SOURCE_PINS.json) make the result reviewable.

```bash
python3.12 -B 04-operator-evolution/verify_r8.py --rkf-root ../Recognition-Kernel-Framework
```

Use the RKF commit listed in the pins, Python 3.11 or 3.12, and Node. These are algebraic proofs and exact computational certificates; a new Lean formalization or physical validation is not claimed. R8 credits the existing compression, observer and KIR theorems rather than duplicating their engine.

## Latest review: the existing EMK master and advanced results

The [EMK master tensor review](02-relational-response/EMK_MASTER_TENSOR_REVIEW.md) corrects the earlier roadmap: master/time, metric/Christoffel/Riemann, holonomy/sheet memory, information tensor and finite cohomology results already exist in the connected repositories. All 191 tests across eight public certificate families pass. Three additional checks identify a deck/compensator sign mismatch, a missing nondegenerate-domain condition, and a determinant-factor interpretation that needs correction. The [review verifier](04-operator-evolution/verify_emk_master_review.py) and [record](04-operator-evolution/EMK_MASTER_REVIEW_VERIFICATION.json) reproduce these findings against unchanged sources. R8 extends observation/curvature contracts while these upstream consistency repairs remain identified and unapplied.

## R7: typed tensors on an admitted EMK carrier

[R7: A typed tensor foundation for EMK](02-relational-response/EMK_TENSOR_CALCULUS_R7.md) develops tensor algebra, transport, derivatives, seam blocks, finite integration, curvature and metric compatibility on the existing native KIR carrier. The [Vault source audit](02-relational-response/EMK_VAULT_SOURCE_AUDIT_R7.md) distinguishes existing work, corrected claims and remaining assumptions. It indexes 143 source files, including all 43 recovered EMK core LaTeX files, with explicit reading depth rather than claiming every file was read in full.

Three load-bearing distinctions are now precise:

- Even-rank tensors cannot distinguish carrier transport `I` from `-I`; a metric can return while a vector does not.
- On the standard real two-mode carrier, both continuous R and K flows cannot preserve a nonzero fixed symmetric metric. Their shared traceless transport preserves an alternating area form; metric selection remains a separate task.
- Tensor return, full carrier return and independent sheet winding are different closure checks.

The fixed-metric obstruction concerns that standard two-mode continuous action; it does not exclude the existing EMK geometry or information metrics. The [exact implementation](04-operator-evolution/emk_tensor_calculus.py) passes 16 focused tests, including all 30 rank/variance signatures at ranks 1–4. Three unchanged public EMK source modules are consumed directly. The verifier also preserves all recorded R1–R6 hashes:

```bash
python3.12 -B 04-operator-evolution/verify_r7.py --publications-root ../Publications
```

See [source pins](04-operator-evolution/R7_SOURCE_PINS.json) and [verification evidence](04-operator-evolution/R7_VERIFICATION.json). Classical tensor/bundle calculus already has rotation and holonomy; R7 specifies the EMK representation and recognition/ledger policy rather than claiming those general operations are new. Follow the [master review](02-relational-response/EMK_MASTER_TENSOR_REVIEW.md) to connect this admitted carrier to existing RTC and metric sectors before choosing further extensions.

## R6 recovery with coefficient error bars

[R6](03-lambda-reference/INTERVAL_PROFILE_RECOVERY_R6.md) encloses recovered depth couplings when the calibrated response coefficients come with supplied uncertainty intervals. It returns the certified prefix, stops at an undecided positive product, and distinguishes that situation from a coefficient box incompatible with any positive profile of the requested depth.

With absolute coefficient errors of `1/100`, R5's original and changed sources still yield disjoint third-cell ranges. Future response bounds now retain both coupling uncertainty and the unknown tail. Their extrema are checked at eight corners using the unchanged native solver.

```bash
python3.12 -B 04-operator-evolution/verify_r6.py --rkf-root ../rkf-r4
```

The 10 focused checks and preserved R1–R5 hashes are recorded in [R6_VERIFICATION.json](04-operator-evolution/R6_VERIFICATION.json). The [implementation](04-operator-evolution/interval_profile_recovery.py) uses exact rational intervals. Observation-error bounds and gain calibration are inputs; obtaining them from an experiment remains a separate task.

## R5 recovery of an unknown depth prefix

[R5: What a boundary response determines about its depth profile](03-lambda-reference/PROFILE_RECOVERY_R5.md) adds an exact inverse procedure without assuming a two-cell period. With a calibrated common product multiplier `z`, the first `m` odd coefficients of the response at `z=0` determine precisely the first `m` positive cell products. For example, `[3,-18,216]` recovers `[3,2,3]`; deeper cells remain unknown.

R5 also proves the limitation at every finite depth: the same weak-response coefficients **and the same exact closed bond** can come from different completed profiles. The explicit sources `(3,2,3,2,...)` and `(3,2,2,1,2,1,...)` both close at `z=1` and share coefficients `3,-18`, but differ at the next coefficient and at a held-out native probe. An additional bound on the unseen tail gives a certified prediction interval.

This applies classical continued-fraction inversion to the existing native source. It extends the inspected fixed-period identification packets by treating arbitrary depth dependence, with a different calibration contract; it does not claim a new general inversion method or a selected physical constant. The [proof's lineage comparison](03-lambda-reference/PROFILE_RECOVERY_R5.md#7-lineage-and-what-is-added) credits the existing RKF and Publications results.

The 14 focused checks pass, including nine comparisons with the unchanged native finite-inverse solver and four source proof replays. Reproduce with Python 3.11 or 3.12 and the same pinned RKF runtime used by R4:

```bash
python3.12 -B 04-operator-evolution/profile_recovery.py '[3, -18, 216]'
python3.12 -B 04-operator-evolution/verify_r5.py --rkf-root ../rkf-r4
```

See the [implementation](04-operator-evolution/profile_recovery.py), [source pins](04-operator-evolution/R5_SOURCE_PINS.json) and [verification record](04-operator-evolution/R5_VERIFICATION.json). Exact coefficients are inputs; recovery from noisy experimental readings is not yet certified.

## R4 native source bridge and closure checks

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

R1 specifies a scalar model, its reference changes, and its complete return invariants. R2 instantiates the existing cut-graded return and compression calculus. R3 supplies an explicit conditional seam/bond realization. R4 connects it to an existing native exact-cut response and checks closure without discarding its amplitude or aperture dependence. R5 determines the recoverable depth prefix from calibrated response coefficients and proves that a finite jet plus exact closure still leaves deeper source freedom. Next supply an independent source law, or an observation protocol with warranted calibration and error bounds. Invariance, identification from data, and universal value selection remain separate claims.

The name `lambda` currently serves several roles. Use `lambda_*` for the fixed reference, `chi_XY` for a response derivative, `gamma` for a decay rate, and `epsilon` for observational resolution until a theorem relates them.

## Provenance

[MANIFEST.json](MANIFEST.json) records the original source checksum and the line ranges of the excerpts. R4 and R5 execute a pinned upstream implementation from a separate checkout; its source files are referenced rather than duplicated here. [CROSS_REPO_LINEAGE.md](CROSS_REPO_LINEAGE.md) records the R1–R4 assessment; the R5 proof adds its own scoped comparison without changing that earlier evidence.
