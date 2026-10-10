# Stopping-rule audit — reuse endpoints, distinguish scope from obstruction

10 October 2026. Owner instruction: inspect whether earlier progress was
stopped by an overly restrictive certificate; reuse valid results and
extend a promising endpoint instead of restarting.

This is a targeted audit, not a replay of every research certificate.
Read sources are pinned in SOURCES.json. No frozen packet was changed.
The old discussion attachments are navigation material, not proof or
authoritative source code.

## Findings from actual sources

| Source / endpoint | What the inspected source actually requires | Classification and action |
|---|---|---|
| Publications YM93 | `theta_value` rejects theta>1/100; the all-coupling gradient helper separately accepts every nonnegative theta. Gap, diagonal-Hessian and continuation constants were proved only in the declared interval. | A valid frozen-domain guard. Failure outside its window means **uncovered**, not zero gap or impossible continuation. Do not remove the guard; prove a successor estimate. The source's radius-two SU(2) normalization must be matched before comparing later YC bounds. |
| YM93 certificate metadata | Model includes `positive_four_dimensional_continuum_mass_gap=False`, while `run()` returns `verdict='PASS'` and the claim status explicitly says continuum open. | No evidence that a missing final classical target prevented this native/lattice packet from passing. The False field records scope; it is not the packet verdict. |
| NSB2 T1 / T2 | Powers of one symmetric 2x2 response remain in span{I,H}; its Hessian-square rigidity theorem assumes a fixed affine calibration, positive Hessian, constant coefficient and identities on an open set. | Genuine obstructions **within those operations and hypotheses**. They do not prohibit new source channels, transported pairings or all nonlinear constructions. A coordinate change alone cannot negate the theorem's correctly transported content. |
| NSB2 spectral reconstruction helper | Positive simple spectrum, calibrated frame and every required distinct-index cubic component are enforced. | A valid information/domain guard for that adapter. At collisions use direct response jets or additional measurements; do not silently divide by zero or declare nonlinear geometry impossible. The frame-free curvature routine is a different adapter. |
| AG1 to TW1 | AG1's turn inequality was initially checked at selected couplings; TW1 provides a subsequent all-coupling proof on that heat-record carrier. | A successful example of retaining the endpoint and strengthening its rule. It still supplies no exact 4D Wilson block law or carrier identification. |
| YC31 to YC32 | Per-source locality did not imply a collective norm. YC32 retained the pulled-back source vectors and used product-vacuum orthogonality / overlap counting. | A prior real repair of an overly coarse **estimate**, not a faulty certificate. Reuse its source-Gram principle for YC37's boundary ports. |
| YC36 | Nearest-pair support distance loses unmatched spectators; Hausdorff distance and both-sided localization give bounded support degree at fixed count. Constants still grow with count. | Geometry of the selected reading was repaired already. A fixed-count theorem is not an all-count or all-word bound. |
| YC37 onsite reduction | The returned metric is retained: the pencil is K-zM+O(z^2/b), not K-zI. | Moving pairing is already implemented. Requiring the raw coefficient matrix to preserve an untransported Euclidean metric would be an additional, unjustified rule. No such erroneous rejection was found here. |
| YC37 join | Exact spectral convolution and ordered boundary histories are retained. A static effective Hamiltonian alone fails even for independent joins. | A real closure requirement. Enlarging the retained record opens algebraic composition, but the extensive Dyson bound does not give uniform iteration. |
| YC37 final boundary warning | A crude squared boundary norm divided by a hidden floor gives O(b^(1/3)). It is explicitly called a diagnostic, not a no-go. | Promising endpoint. YC38 improves the same selected centered vacuum-source estimate to O(b^(-1/3)) using inherited locality and the full Gram. The original all-input/all-history gate remains. |

**Diagnosis:** no incorrect rejection was established in the inspected
executable guards. A recurrent risk is interpreting a sufficient scalar
bound or a selected representation as the boundary of all admissible
constructions. The latest stages already correct several such losses.
There is also concrete room to improve an estimate by retaining old
source information: YC38 supplies one such extension.

## Rules applied in this audit

1. State the carrier, physical pairing, source class, allowed operations,
   domain, energy units and target before testing a claim.
2. Record **established**, **uncovered by this rule**, or **obstructed in
   this specified class** separately. A failed sufficient bound is not a
   mathematical no-go, and a missing continuum bridge does not erase an
   established lattice theorem.
3. Carry the entire source Gram and returned metric before replacing
   them by a scalar norm. Preserve cross terms, signs and hidden returns.
4. Under a frame change transport every reference, cut, source and pairing.
   Under a physical change derive an evolution or comparison law instead
   of assuming the old invariant value stays constant.
5. For iteration record which data close under the next operation.
   Two-point response, all energy derivatives, and all ordered boundary
   words are different contracts. Their quantifiers cannot be exchanged.
6. Reuse frozen results as inputs. A new family or wider domain receives
   a successor proof with new hypotheses; do not weaken the old guard or
   rewrite historical evidence to obtain a PASS.
7. Follow the inherited theory-first preference: written proof and source
   audit now, additional executable/numerical certification later. Label
   actual verification honestly. These rules are research bookkeeping,
   not a requirement to pass the final target before exploring.

## What changed and what benefit was obtained

The new note joins two existing endpoints: the size-independent locality
of the vacuum frame (YC31-32) and the cb hidden floor of a spatial count
cut (YC37). It constructs the boundary-port Gram on the actual block
vacuum, proves its summable row bound, and then uses operator order after
the potentially nonlocal count projection. Projection is not assumed
to preserve entrywise covariance decay.

For N_boundary local ports and bounded scalar coefficients, the old
triangle bound is g^2 N_boundary^2/(cb). The new bound is
g^2 N_boundary min{N_boundary,C_Gamma}/(cb). For regular 3D blocks this
has a decaying O(b^(-1/3)) branch, and its lifted-response norm squared
has an O(b^(-4/3)) branch. The practical crossover constant is not
evaluated; the min form prevents claiming an unproved finite-size gain.

This is a **source-scoped theoretical improvement**, not a new measured
mass gap, a larger coupling window or a complete spatial contraction.
The exact product-vacuum Z-port example in YC38 shows why a zero source
Gram can coexist with a large full operator. Thus the remaining task
cannot be removed just by changing certificate labels.

## Reading / verification scope

- Read YM93's theorem and actual guard, model metadata and continuation
  helper code; read NSB2's theorem and spectral reconstruction guards.
- Read the existing AG1/TW1 construction and its stated carrier boundary;
  read YC31-32, YC36-37 and the preceding target-reference audit.
- Used the current local YC37 branch, newer than the uploaded discussion's
  YC36 endpoint. No claim to have reread every YM1-93 or RH proof.
- Derived YC38 using the stated inherited locality and gap inputs; checked
  its covariance expansion, Schur order, source transformations, unit
  scaling and zero-reference cases mathematically.
- Did not run old certificate suites or create a new numerical/theorem
  certification packet. No claim of formal proof, independent review,
  a computed C_Gamma or all-history closure.

Next: an operator-valued connected boundary-history bound with controlled
dependence on input/source count and word length. This is narrower and
more informative than simply requesting another fixed-circle PASS.
