# When does a reference lambda become a determined constant?

This note contains elementary consequences of the displayed formulas in the supplied draft. It does not add an empirical result or claim to derive a physical constant.

## 1. Fixed during evolution is not uniquely selected

A model can hold `lambda_*` constant while its states and response coefficients vary. That statement alone does not establish which constant must be used or whether two different references describe distinguishable systems.

## 2. Exact rescaling freedom for unrestricted scalar coefficients

**Proposition.** Let `lambda_*>0`, let each `k_i` be an integer, and suppose that the coefficients `n_i` are unrestricted real or complex scalars. If the model's observations and constraints depend only on

`X_i=n_i lambda_*^{k_i}`,

then those observations and constraints do not determine `lambda_*` uniquely.

**Proof.** For any `c>0`, set `lambda_*'=c lambda_*` and `n_i'=c^{-k_i}n_i`. Then

`n_i' (lambda_*')^{k_i} = c^{-k_i}n_i c^{k_i}lambda_*^{k_i} = X_i`.

Thus every permitted `c` gives the same `X_i`. The assumptions ensure that the adjusted coefficients are still allowed and that no other condition distinguishes the references. For `c != 1`, the reference changes. QED.

**Boundary of the result.** It does not apply unchanged when coefficients have a separately fixed lattice, integrality condition, or native constraints that are not preserved by rescaling. Those conditions must be supplied and their consequences proved. Nor does it exclude invariant dimensionless combinations of quantities or independent calibration.

## 3. A uniformizer has its own reference freedom

Let `F` already be a discretely valued field, let `v(pi)=1`, and set `U={u in F: v(u)=0}`. Every nonzero `x` has a unique decomposition `x=pi^k u` relative to the chosen `pi`.

For any `a in U`, the replacement `pi'=a pi` still has valuation one, and

`x=(pi')^k (a^{-k}u)`.

The new coefficient is still in `U`. Therefore the valuation and decomposition alone do not distinguish `pi` from `a pi`. Additional structure could distinguish one; it is absent from the decomposition itself. This follows from the standard uniformizer framework described in [the Stacks Project](https://stacks.math.columbia.edu/tag/00P7).

## 4. The shift-scale operator does not select a numerical lambda

Use the explicit model on `ell^2(N; C)` with `S e_i=i e_i` and `L e_i=e_{i+1}`. On finitely supported vectors, direct evaluation gives `[S,L]=L`.

For any positive scalar `a`, define `S_a=aS`. Then

`[S_a,L]=aL`.

Every `a` produces such a ladder. Rewriting its spacing as `lambda_*` therefore does not derive a distinguished spacing. The model would need a native condition that selects one equivalence class or an invariant ratio, with no adjustable parameter reintroduced elsewhere.

## 5. A concrete research target

The next useful theorem would specify:

1. A native algebra and allowed changes of reference or representation.
2. A precisely defined dimensionless quantity `r` invariant under those changes.
3. An equation or structural characterization for `r` derived from the native relations, with no freely tuned coefficients used to set the answer.
4. Existence and uniqueness in a stated admissible class.
5. A map from `r` to a response observable that distinguishes alternative models or values.

A ratio of nonzero operator eigenvalues, under a common scale change, illustrates the kind of quantity that survives simple rescaling. A particular native operator, its spectrum, and its physical interpretation would still have to be constructed. This is a candidate direction, not an established derivation in `4ways.tex`.
