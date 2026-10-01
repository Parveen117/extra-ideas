# Way-2: Relational and response numbers

Source: the section of that name in [4ways.tex](../4ways.tex). [source-excerpt.tex](source-excerpt.tex) preserves it verbatim as a LaTeX fragment.

## Core idea

Describe a system through responses and relations rather than only its coordinates. Write `chi_XY=dX/dY` for a response and reserve `lambda_*` for the proposed fixed reference.

## What can be used now

On an admissible differentiable curve with `dY/ds != 0`,

`chi_XY = (dX/ds)/(dY/ds)`.

This ratio is invariant under a regular reparameterization of the same curve. That useful invariance does not make the ratio independent of the state or of the curve. In a multivariable system, state which variables or constraints are held fixed.

## What needs development

Specify the admissible transformations and show when `Gamma=chi_p/chi_v` is invariant. Requiring that an evolution preserve a ratio is a constraint; deriving that preservation is a separate result. State the domain where denominators are nonzero and provide a justified treatment of any seam endpoint.

For a proposed Jacobian, check the compatibility conditions needed for it to arise from actual observables. Do not equate products of instantaneous derivative responses with a finite sequence of state transformations without an integration argument.

## Next target

Continue from the existing EMK master before physical response laws. The [master review](EMK_MASTER_TENSOR_REVIEW.md) credits the already developed tensor/time, metric/curvature, holonomy, information and topology layers, and identifies three specific consistency corrections. [R7](EMK_TENSOR_CALCULUS_R7.md) supplies typed E/M pairing, tensor increments, ordered integration and curvature on one admitted carrier. [R8](CURVATURE_OBSERVATION_R8.md) certifies full curvature, compressed readout and recomputed reduced curvature. [R9](CURVATURE_BALANCE_R9.md) now supplies a conditional native balance mechanism with an exact information/branch ledger. Next identify its cut grades, branch record and recognition-event weights in a chosen RTC sector while retaining derivative compatibility and the existing sheet policy.

## R9 conditional information and curvature balance

[CURVATURE_BALANCE_R9.md](CURVATURE_BALANCE_R9.md) derives the exchange-balanced endpoint from the existing cut-swap, proves which curvature grades cancel or survive, and gives an exact event-by-event imbalance law. It transforms the complete information family and its derivatives, records the visible information loss and exact branchwise recovery, and distinguishes recognition curvature from Riemann curvature. [R9_VERIFICATION.json](../04-operator-evolution/R9_VERIFICATION.json) binds 24 exact tests and nine symbolic native replays. Off-diagonal entries and observation in general are not treated as universal curvature or flatness criteria.

## R8 curvature and response development

[CURVATURE_OBSERVATION_R8.md](CURVATURE_OBSERVATION_R8.md) proves the exact discarded-sector correction for observation, constructs counterexamples in both directions, and recovers every current KIR odd-curvature target from two fixed cross-sector channels. Dimension restriction follows an explicit two-form pullback. Onsager symmetry is equivalent to flatness for the stated constant linear scalar response; a variable symmetric response can still be curved. The proof keeps active cut sign reversal, passive reference changes and independent sheet memory distinct. See the [certificate](../04-operator-evolution/R8_VERIFICATION.json) for 23 exact tests and native proof replay.

Role in the four-route program: connect an algebraic or operator result to a quantity that can be compared across states or measured. See [the shared assessment](../ASSESSMENT.md).

## R1 development

[RESPONSE_RETURNS_R1.md](RESPONSE_RETURNS_R1.md) identifies closed-return products and equal-endpoint path ratios that are independent of local references. It also explains why ordinary derivatives along one common tangent give identity returns and why conservation is an additional condition.

## R2 development

[PROJECTED_RESPONSE_R2.md](PROJECTED_RESPONSE_R2.md) shows how a fixed full-return sign structure can coexist with a variable response after a bond projection. The exact defect is algebraic; a nonnegative norm interpretation requires the additional metric assumptions stated there.

## R3 development

[TYPED_SEAM_R3.md](TYPED_SEAM_R3.md) separates the scalar constraint from a state operator, then fixes the projection by an explicit rule for its removed subspace. It gives the compatible return response and explains why scalar closure alone can hide leakage.

## R4 source connection

[R4](../03-lambda-reference/NATIVE_BOND_BRIDGE_R4.md) consumes the existing paired-depth response and carries its raw amplitude and finite-tail enclosure through the bond test. It demonstrates with actual source responses why normalization and one finite closed aperture cannot certify completed closure.

## R5 inverse response and remaining uncertainty

[R5](../03-lambda-reference/PROFILE_RECOVERY_R5.md) reconstructs the first `m` positive cell products from `m` calibrated odd weak-response coefficients. It proves the first order at which a deeper change becomes visible and uses the unchanged source's transfer map to enclose future responses with the unknown tail retained. Unknown input and output gains have an explicit effect on the recovered products; noisy coefficient recovery remains a separate task.

## R6 response error propagation

[R6](../03-lambda-reference/INTERVAL_PROFILE_RECOVERY_R6.md) carries supplied coefficient error bars into depth-coupling intervals and a future response interval. It retains the unobserved tail and uses the native solver to check both prediction extrema. Calibration and observation-error validity remain part of the input contract.
