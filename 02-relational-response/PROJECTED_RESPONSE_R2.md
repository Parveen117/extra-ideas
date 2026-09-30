# R2: A fixed full return can have a variable retained response

The [R2 proof](../03-lambda-reference/AGHORA_RETURN_R2.md) specializes the source's exponential flow to its Aghora-odd generator and constructs `R_t=A exp(tG)`. The full return obeys `R_t^2=I`, with eigenvalues restricted to `+1` and `-1`.

A retained observation is represented here by a specified projection `B`. The corresponding response is the compression `C=BR_tB` on `ran B`. It is not automatically an involution. Its exact discrepancy is

`B-C^2 = BR_t(I-B)R_tB`.

This compares two full returns with two returns that include an intermediate projection. The omitted component can contribute to the second full return. An intervening cut prevents that contribution from returning.

If `B` is orthogonal and `R_t` is self-adjoint, the discrepancy is a nonnegative squared-norm operator. In the exact example documented in R2, the full return sends `(1,0)` to `(4/5,3/5)`. Projection onto the first coordinate leaves `(4/5,0)` and removes a squared norm of `9/25`. Two full returns restore `(1,0)`; a projection after each pass leaves `(16/25,0)`.

These are conditional consequences of an explicit chosen model. The fractions are not universal constants and the norm balance is not an entropy definition. Selecting a particular response requires deriving the retained subspace or projection from the native rules.

A constant full-return sign and a variable projected response can therefore coexist. The bond-selection law is the next unresolved connection between the source's seam terminology and a unique observable value.
