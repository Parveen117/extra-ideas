# R15 — Cut, recognition, native iota and topology

Research owner: Monty Dabas. Corrective development: 1 October 2026.

**The foundation proceeds from cut records and operators. A complex scalar is a later chart.** R14's multiplication by `a+iota b` was a scalar specialization and was incorrectly promoted as the foundation of the full UGD/recognition program. Its calculations do not determine the scope or limitations of the full native theory. See [the explicit R14 correction](R14_SCOPE_CORRECTION.md).

This packet audits all **48 theorem/lemma/proposition/corollary environments** in the supplied `emk_topology.tex`, repairs its recognition/topology construction, and constructs a native quarter-turn from signed cyclic cut records. It retains scale and seam information. It does not certify the upload's claim that every later mathematical and physical structure follows from its three prose axioms.

The original is preserved [byte for byte](../04-operator-evolution/certificates/r15/emk_topology.original.tex), SHA-256 `55f449e7eb3a2bde9e77193f595a6a5eeff3923da2bea076d00f0e9907b8cebd`. The [claim ledger](../04-operator-evolution/certificates/r15/CLAIM_LEDGER.json) identifies the exact source environments, judgments and replacement results. The [corrective LaTeX edition](emk_topology_verified_R15.tex) contains the central proofs. “Verified” here means scoped written proofs and exact computational checks, not a proof-assistant certificate or approval of every original claim.

## 1. What belongs at the source

The upload names Śūnya, a first cut `kappa`, and a self-cut `chi`. Its first-cut arrow has domain Śūnya; the later recursion applies the same symbol to arbitrary appearances. Its self-cut initially has a trace domain, then also acts on every appearance. These extensions require a type declaration. Composition cannot silently repair a domain mismatch.

R15 uses typed finite cut records first. When a common appearance carrier is used below, its two total actions and the self-cut output record are explicitly part of that presentation. This is a precise interpretation of the manuscript's recursion, not a proof that the original prose uniquely selects this presentation. An arbitrary set of labels, complex coordinates, topology, metric, stochastic law, physical clock, and thermodynamic potential are not primitive inputs to the construction.

The proof language uses finite words, records, functions, subsets, logical equality and induction. These are **metamathematical infrastructure**, openly used also by the upload's own set notation and recursion. Object-level recognition equality can be constructed; logical equality cannot be eliminated merely by renaming it recognition. No classical analytic theorem is imported into §§2–8. This does not mean that standard logic/set/function reasoning has been derived from the three prose axioms.

## 2. R15.1 — Typed paths and the free source test

Give every cut occurrence its source and target. An empty path exists at each object. Paths compose by joining matching endpoints; concatenation is associative, and the empty paths are the separate local identities. Evaluation into any other system of maps with those same signatures is uniquely determined by the images of the generating cuts.

**Proof.** A length-zero path must evaluate to the corresponding local identity. A path of positive length is its shorter prefix followed by its last arrow, so its evaluation is forced recursively. Associativity follows because both parenthesizations produce the same ordered list. This proves existence and uniqueness by induction. No coordinates or numerical field occur.

For the manuscript's common-carrier interpretation, finite words in `kappa,chi` applied to `p0` give a free term construction. No relation equates the nonempty word `chi kappa chi kappa` to the empty word. Thus merely composing the two acts does **not** prove `(chi kappa)^2=I`. The unchanged RKF symbolic engine independently replays this distinctness in the free presentation.

There is also a three-state witness: extend the appearance cut as `kappa=(1,2,0)` and use `chi=(0,1,2)`. Recognition by exact self-cut output distinguishes the states. Then `J=chi kappa` sends `0→1→2`, so `J²(0)=2`, not `0`. A separate initial arrow may still send the ground to the seed `0`. This tests the manuscript's algebraic totalization; it is not a theorem about the metaphysical interpretation of Śūnya.

Consequently a four-cycle or inverse-cut law must be obtained from a specified native construction. It cannot be inserted under the word “derived.” This finding does not refute an existing, stronger native iota construction elsewhere in the repositories.

## 3. R15.2 — Recognition that survives every future cut

Let `G` be the cut actions on the explicitly typed appearance carrier `X`, and let `r(x)` be the recorded self-cut output. For the exact-output interpretation of the upload, take `r=chi`. If a coarser readout is used, record that readout explicitly. Define

\[
x\approx y\quad\Longleftrightarrow\quad
r(w x)=r(w y)\quad\text{for every finite admissible word }w,
\]

including the empty word. The common-carrier statement below uses total actions. For partial typed actions, equal source type and equal future availability must also be retained; the executable exhaustive check covers the total case only.

**Theorem.** The relation `≈` is an equivalence, is preserved by every cut, and is the greatest cut-compatible equivalence contained in current recognition `ker(r)`.

**Proof.** Reflexivity, symmetry and transitivity hold separately for each equality of records, hence for their conjunction. If `x≈y`, then for any generator `g` and future word `w`, the word `wg` is also a future, so `r(wgx)=r(wgy)`. Thus `gx≈gy`. Finally, any cut-compatible equivalence contained in `ker(r)` preserves every finite word by induction, and therefore implies `≈`.

This is the most information-losing quotient that is still faithful to all these future records. It is not full history equality: a discarded history component must remain in the observation family if future equality is meant to retain it. Adding seam, winding and scale records refines the equivalence automatically. On the completely free term carrier with literal self-cut outputs, no distinct terms collapse; nontrivial recognition requires actual source identifications or a declared readout.

For a finite carrier, begin with `R0=ker(r)` and repeatedly retain only pairs whose images under every generator remain related. Each strict refinement removes a pair; the process terminates. Its stable relation equals `≈` by the theorem. A second implementation enumerates the entire finite transition monoid and intersects the corresponding record kernels. R15 compares these independent algorithms exhaustively.

## 4. R15.3 — A correctly typed Eye and seam

The Eye is the quotient map `q:X→Q=X/≈`. Its square is not an endomap expression: `q(q(x))` has the wrong domain. Instead define saturation on appearance predicates,

\[
\operatorname{Sat}(U)=q^{-1}(q(U)).
\]

**Theorem.** Saturation is extensive, monotone, union preserving and idempotent.

**Proof.** It includes exactly the equivalence classes meeting `U`. A class already included is unchanged by another saturation. Inclusion and union statements follow from the same membership test.

Cuts descend by `ḡ([x])=[gx]` because of R15.2. Define the corrected zero-change seam

\[
S_g^0=\{x:q(gx)=q(x)\}=q^{-1}(\operatorname{Fix}(\bar g)).
\]

This is a fixed-point statement on the quotient with consistent types. It differs from the upload's seam, which compares `chi(Cx)` with `chi(x)` **again under recognition**. For `R=ker(chi)`, the upload's seam tests `chi²(Cx)=chi²(x)`, whereas its zero healing cost tests `chi(Cx)=chi(x)`.

An exact witness has `chi=(0,0,1)`, `C=(2,2,2)` and `x=0`. The original seam test passes but healing cost is **one**, not zero. More generally, for this exact-output interpretation, equality persists after further `chi`, so the original seam is the set of costs at most one. R15 does not silently redefine the original cost to manufacture a pass.

## 5. R15.4 — Topology derived from saturated cut stability

Define

\[
\tau_{\mathrm{cut}}=\{U\subseteq X:
\operatorname{Sat}(U)=U,\quad g(U)\subseteq U\ \forall g\in G\}.
\]

**Theorem.** This is a topology, in fact closed under arbitrary intersections. It is exactly the inverse image of the forward-invariant topology on `Q` and exactly the family of upper sets of the preorder generated by recognition steps and directed cut steps.

**Proof.** Empty and full sets satisfy both conditions. A union of saturated sets is saturated, and a cut image of a union lies in that union. For an arbitrary intersection, saturation follows because all memberships are constant on equivalence classes. If `x` lies in every member, each `gx` lies in every member as well. Thus all intersections are stable. Saturated sets are exactly inverse images under `q`; cut stability translates to stability under the descended actions. Finally stability under single recognition/cut steps is equivalent, by finite induction, to stability under their reflexive transitive closure, which is the stated preorder.

This proves a concrete cut topology without choosing a manifold or Euclidean metric. It does **not** assert that every arbitrary generator is continuous in this topology when generators do not commute. It is also not the “coarsest topology making endomaps continuous”: the indiscrete topology already makes all endomaps continuous. An initial topology needs specified observable maps and codomain topologies.

The upload's weaker condition `q(C(U))⊆q(U)` does not force saturation or intersection closure. With `chi=(0,0,2)`, `kappa=(2,2,0)` and `R=ker(chi)`, the sets `{0,2}` and `{1,2}` satisfy the condition for **every generated word**, but their intersection `{2}` does not. The upload's separate sufficient assumption that the quotient's direct image preserve all finite intersections would fix closure, but is not derived. Indeed such preservation on all subsets forces the quotient map to be injective: two different points in one fibre violate it using their singleton sets. The repair permits nontrivial recognition fibres.

## 6. R15.5 — Counting cuts, distance and event records

The cut count of an empty record is zero; adjoining one occurrence forms its successor. Concatenation adds these counts by induction. This derives an event index and word length from the record syntax; it does not supply physical time or a real-valued clock.

On `Q`, let `d⁺(a,b)` be the least number of cut occurrences in a word carrying `a` to `b`, with infinity when no word does so. Then

\[
d^+(a,b)=0\iff a=b,\qquad
d^+(a,c)\le d^+(a,b)+d^+(b,c).
\]

**Proof.** A zero-length word is the identity. Two finite realizing paths concatenate to a path of summed length, giving the triangle inequality; if a right-hand term is infinite the extended inequality is automatic. A nonempty subset of finite cut counts has a least member by induction on counts.

This is a directed counting distance, not necessarily a symmetric metric. The chain `0→1→2→2` has forward distance two from `0` to `2` and no return. It is also a new, explicitly identified construction, not the upload's healing cost. A scalar metric, stochastic noise law, continuous angle and physical event actualization have not been smuggled into it.

## 7. R15.6 — Native iota from the UGD phase operator

The existing UGD-1 source already carries phase as a cyclic exponent, alongside scale and a seam ledger. Reuse that native phase alphabet rather than starting with complex numbers. Its size is a declared native parameter; A0–A2 in this upload do not select it.

Take a phase alphabet of size `4m`, where `m` is a positive cut count. Form the free signed cut ledger on phase marks `e_j`, indices modulo `4m`. “Signed” means formal records modulo cancellation of opposite occurrences; integer coefficients are counts in this construction. Let the phase successor be `U e_j=e_(j+1)`. Put

\[
f_j=e_j-e_{j+2m}\quad(0\le j<2m),\qquad
V^-=(I-U^{2m})G_{4m},\qquad
\iota=U^m|_{V^-}.
\]

**Theorem.** The `f_j` form a free signed basis of `V⁻`; on this sector, `iota²=-I` and `iota⁴=I`. Define reflection `K e_j=e_(m-j)`. It preserves `V⁻` and obeys `K²=I` and `K iota=-iota K` there.

**Proof.** Every image `(I-U^(2m))e_j` is one of `f_0,…,f_(2m-1)` or its negative, proving spanning. In a zero combination, its coefficient on each first-half mark is precisely the corresponding coefficient, proving independence. A half-cycle exchanges the two terms of every `f_j`, so `U^(2m)f_j=-f_j`; hence the two power identities. Reflection twice returns every phase mark. It reverses successor: `K U^m=U^(-m) K`. On `V⁻`, `U^(-m)=-U^m`, giving anticommutation. These are identities of signed cut records, before any matrix or complex chart.

At `m=1`, `iota f_0=f_1`, `iota f_1=-f_0`, `Kf_0=f_1`, `Kf_1=f_0`. The familiar two-component action and the expression `a+iota b` can be written **after** this operator construction. Their scalar chart is not the source of the generalized numeral.

This construction preserves its scope: the whole phase module is larger than `V⁻`, and `U^(2m)e_0=e_(2m)≠-e_0` there. Scale indices and seam ledger components remain separately carried. Extending `iota` by the identity on a seam ledger would square to `(-I,I)`, not minus the identity on that enlarged object. Thus this proof does not flatten UGD into a complex scalar or assert a global scalar action on every linked-seam structure.

The construction explains exactly where the quarter-turn law comes from: the cyclic phase action and its signed half-cycle sector. It does not claim that naming first-cut and self-cut alone proves cyclicity, selects `4m`, or identifies either act with `K` or `iota`. Those identifications need the source's typed operator lineage. No new classical analytic premise is used here.

## 8. R15.7 — Closed endpoints retain path memory

In the four-phase presentation, a positive path has count `n=4w+r`, `0≤r<4`, by repeated removal of complete four-step blocks. Its visible phase is `r` and its complete-loop ledger is `w`. Concatenation obeys

\[
(r,w)\cdot(s,v)=
\bigl((r+s)\bmod4,\ w+v+\lfloor(r+s)/4\rfloor\bigr).
\]

**Proof.** Add the two counts, remove all full blocks, and keep the remainder. Associativity follows from concatenation/count addition, not from erasing the carry. A zero path and a four-step path have the same endpoint but loop ledgers zero and one. The latter is not an empty path in the history construction.

Consequently the operator identity `iota⁴=I` is compatible with nontrivial path memory: operator evaluation forgets information. A source equation for logarithmic increments can be path independent on the lifted record while differing on two paths with the same projected endpoint. The upload's global same-endpoint argument erased this distinction. No classical contour integral is needed to exhibit the lost information.

## 9. R15.8 — Involution, compass and geometry are different objects

For any `J²=I`, induction gives `J^(2n)x=x` and `J^(2n+1)x=Jx`. Thus its endpoint orbit has at most two members. It cannot by itself produce a four-point orbit, an unbounded layer tower, or an infinite helix. If also `J²=J`, then `J=I`. An involution is therefore not generally an idempotent Eye. A two-state swap has no fixed point at all.

The upload's chosen complex chart `J(z)=i conjugate(z)` is the real-linear swap `(x,y)↦(y,x)`. Its fixed set is the whole diagonal, including `(-1,-1)`, not just the positive ray. This checks the upload's explicit representation; that representation is not an input to the topology proof. Also `[J,J]=0` in any associative operator algebra. Nonzero mixed cut-order curvature can instead involve two distinct operators; the native `K,iota` relation gives `[iota,K]=2 iota K`.

A declared four-edge directed compass has only empty/full forward-invariant subsets, because any included vertex reaches all vertices. On its four distinct vertices this is an indiscrete topology. Collapsing the entire orbit produces a **one-point quotient**, a different carrier. Retaining paths produces a category in which a four-edge loop is not automatically the empty path. None of these facts derives Maxwell identities, potential functions or constitutive thermodynamics. Those require actual response/differential structure, absent from A0–A2 here.

## 10. Repository lineage, with imports stated in each citation

All paths and hashes are pinned in [SOURCE_PINS.json](../04-operator-evolution/certificates/r15/SOURCE_PINS.json). References below credit existing work; they are not blanket certifications of entire repositories.

| Reference and exact scope | What can be reused | Native inputs / imports and evidence |
|---|---|---|
| [RKF native primitive order](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/NATIVE_PRIMITIVE_ORDER.md) | Cut/refinement before signed scalars, quarter-turn, residue paths and later charts | Dependency protocol, not a proof of A0–A2 sufficiency. Native completion and orientation contracts must travel with the result. |
| [RKF N01–N03](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md) | Typed path algebra, separate history/sheet grading, future observer completion | Declared native cut field and action/observer inputs; no primitive classical Hilbert space. N03 is the existing linear analogue of R15.2, not a newly discovered principle here. |
| [RKF automatic regular action WB1–WB4](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/AUTOMATIC_REGULAR_MODEL.md) | Normal basis and faithful operator action generated from relations | Admitted native presentation, termination and overlap audit. Matrices are outputs. R15 replays ten identities/distinctions and the four-element KIR normal basis through this unchanged engine. |
| [UGD-1 finite numeral source](https://github.com/Parveen117/Publications/blob/e1dc4e3773f063f14e56cf30222c8e8a504519cd/papers/emk-ugd-algebra/certificates/ugd1_numerals.py) | Cyclic phase exponents, scale convolution, conserved seam carries, projection blindness | Finite phase alphabet, exact rational scale chart and stated seam rules. No evaluated complex root or transcendental arithmetic. R15 replays the exact original certificate; infinite numerals and general linked-seam composition are explicitly outside that source's certificate. |
| [EMK-T1](https://github.com/Parveen117/Publications/blob/e1dc4e3773f063f14e56cf30222c8e8a504519cd/papers/emk-ugd-algebra/certificates/emkt1_master_tensor_and_time.py) | Derived cut-swap `J=K` in its specified presentation; independent temporal residue | Finite native KIR algebra represented with exact rational arrays; declared active masks and ledger policy. This is not identification of arbitrary upload `chi kappa` with `K`. Cited, not newly replayed in R15. |
| [LAM-2 arithmetic quarter-turn](https://github.com/Parveen117/Publications/blob/e1dc4e3773f063f14e56cf30222c8e8a504519cd/papers/lambda-seam-calibration/README.md) | Exact `tau(chi_4)^2=-4` in a declared cyclotomic residue quotient | Arithmetic presentation and primitive character are inputs. **Classical analytic imports:** the separate theta-seam comparison uses named theta/Poisson anchors. Only the finite arithmetic instance is relevant to iota here; no entire-capsule “pure foundation” claim. |
| [RKF F00/E](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorems/foundation/F00E_NATIVE_EULER_FROM_IOTA_COMPLEX.md) and [F00G](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorems/foundation/F00G_NATIVE_LOGARITHM_AND_POWERS.md) | Native exponential, Euler identities and logarithm after the cut-scalar/iota construction | Completed ordered native scalar field is an explicit input; no imported classical exponential/logarithm is needed for the native theorem. These start **after** the missing source step, so do not prove A0–A2 selects an oriented cyclic operator. Replay evidence remains in R14's pinned record. |
| [RH F00J/F00K branch](https://github.com/Parveen117/RH-Framework/tree/c588ded973a395b5fce37c82f670616160979a2e/theorems/foundation) | Native period, scalar polar form and logarithm branches | Native ordered completion and prior Euler functions. Written results; no new standalone formal/computational certificate supplied by R15. Applicable later to their scalar sector, not a substitute for UGD state data. |
| [Morphic algebra certification status](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/sources/rkf_reference/theorum/morphic_algebra/CERTIFICATION_STATUS.md) | Existing source's own claim boundaries | Explicitly `RNKE_VERIFIED_WITH_OPEN_OBLIGATIONS`, full manuscript certification false: 9 algebraic claims, 10 conditional claims, 37 incomplete claims, one meta guard. A repository “certified” label must not erase these distinctions. |

The R14 cross-repository ledger remains available for RH/Yang–Mills dependencies. A certificate for a scoped lemma is not a certificate for RH, a continuum mass gap, or all mathematics in a manuscript. Nor does an open claim in this upload invalidate a proved result in another native presentation.

## 11. Verification and the remaining exact task

Run from `extra-ideas`:

```bash
python3.12 -B 04-operator-evolution/verify_r15.py \
  --rkf-root /path/to/Recognition-Kernel-Framework \
  --publications-root /path/to/Publications
```

The [verification record](../04-operator-evolution/R15_VERIFICATION.json) checks original-source integrity, complete 48-claim coverage, all pinned dependencies, all 3,678 two-action systems/observation partitions on one to three states, 29,290 saturation cases, nine concrete boundary witnesses, seven rejected mathematical replacements, ten native symbolic replays, a rejected altered native proof, and signed phase-operator bases for six alphabet sizes. The existing UGD-1 certificate is regenerated in memory and compared to its original bytes/hash; no upstream pin is overwritten. R1–R14 theorem/code/certificate evidence remains unchanged.

The **whole-source verdict is NOT CERTIFIED AS WRITTEN**. The corrected constructive results have proofs and scoped checks. The outstanding root obligation is precise: supply a typed native derivation identifying this manuscript's first-cut/self-cut acts with the cyclic/oriented operator construction and its full active records. Cyclicity, closure, inverse existence and analytic completion cannot be inferred from the names of those acts. Where an existing repository result supplies a step, cite that theorem with its actual premises; where it does not, deriving that step remains the task. R15 supplies the operator-before-chart route and stops the scalar example from being mistaken for that missing proof.
