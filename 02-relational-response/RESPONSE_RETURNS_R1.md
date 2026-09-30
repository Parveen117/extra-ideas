# R1: Responses that produce a reference-independent return

The full hypotheses and proofs are in [RETURN_INVARIANTS_R1.md](../03-lambda-reference/RETURN_INVARIANTS_R1.md).

For a process from the one-dimensional space at `i` to the space at `j`, write its nonzero response coefficient as `a_ij`. Under independent coordinate changes `x_i'=g_i x_i`, its value becomes `a_ij'=g_j a_ij/g_i`.

An open-path coefficient depends on the endpoint references. A closed-path product, or a ratio between two processes with the same endpoints, cancels those factors. This is the first explicit response quantity in this development whose value does not depend on the chosen local units.

## A distinction that prevents a false derivation

If all responses are ordinary ratios `(dX_j/ds)/(dX_i/ds)` at the same state and along the same parameter, every closed product is one. A nontrivial return therefore needs path-dependent operations or additional structure. The proposed finite transport network supplies such a structure as an explicit model assumption.

Do not impose all transport equations on one fixed nonzero global assignment: nontrivial returns change an input after the sequence of operations. They are allowed precisely because the sequence is a process, not a static system of simultaneous equalities.

## Constancy test

For positive time-varying coefficients, a return is constant exactly when the oriented sum of the logarithmic rates along its cycle vanishes. All returns are constant exactly when those rates are vertex differences. This distinguishes a mere change of local reference from a change of the actual return data.

## Immediate target for the source chain

The manuscript's six-stage response chain has no independent cycle in this model. Identify and derive a return operation from the intended native rules before treating a cycle multiplier as a selected lambda. An independently adjustable closing coefficient would leave its numerical value free.
