# Way-1: Intrinsic multidimensional numbers

Source: the section of that name in [4ways.tex](../4ways.tex). [source-excerpt.tex](source-excerpt.tex) preserves it verbatim as a LaTeX fragment.

## R14 foundation route

[R14 §§1–3](../02-relational-response/NATIVE_SOURCE_FOUNDATION_R14.md) uses established native cut arithmetic, derives the radial/complement split from dagger, and reads theta from the native source using the existing period/polar/logarithm theorems. The [source ledger](../02-relational-response/DERIVATION_SOURCE_LEDGER_R14.md) distinguishes native constructions from ordinary-field presentations and imported adapters. This route does not repair every claim in the original Way-1 manuscript below.

## Core idea

Represent a state by a vector in `F^n`, where `F` is an explicitly defined field. Use componentwise addition and multiplication. The coordinate idempotents `e_i` satisfy `e_i odot e_j=0` for `i != j` and decompose the algebra into its coordinate factors.

## What can be used now

Once the field is specified, this is a concrete commutative product algebra. It supports coordinate projections, matrix actions, and a precise starting state space. For `n>1` it has zero divisors; it is not itself a field.

## What needs development

The source's underlying field construction is incomplete, and its exact-shell grading fails. Begin with a defined discretely valued field and the filtration

`F^k(F^n) = {x : min_i v(x_i) >= k}`.

Construct the associated graded algebra through quotient spaces rather than declaring exact valuation shells to be graded pieces. Specify which changes of basis are admissible: arbitrary invertible linear maps need not preserve coordinatewise multiplication.

## Next target

Use [R7's admitted carrier and scalar layer](../02-relational-response/EMK_TENSOR_CALCULUS_R7.md#2-admitted-scalar-and-carrier-layer) to keep intrinsic states, scalar coefficients and KIR operators distinct. Tensor slots are now explicitly typed; nonlinear/partial UGD or morphic states still require a representation theorem before they inherit multilinear operations. The next algebraic obligation is that representation and its admissible morphisms.

Role in the four-route program: a clear representation of states and components. See [the shared assessment](../ASSESSMENT.md) for explicit counterexamples and repairs.

## R1 development

[COORDINATE_TRANSPORT_R1.md](COORDINATE_TRANSPORT_R1.md) supplies a finite labelled coordinate model whose edge operators reproduce response paths and closed returns. It uses a specified field; the source's valued-field construction remains a separate problem.

## R2 development

[AGHORA_MODES_R2.md](AGHORA_MODES_R2.md) derives the positive/negative sector split from the source involution. A nonzero Aghora-odd generator needs both sectors, giving a precise reason to extend the scalar model to multiple modes.

## R3 development

[SEAM_COMPLEX_COORDINATES_R3.md](SEAM_COMPLEX_COORDINATES_R3.md) constructs a complex action on two real seam coordinates under an explicit bond rule. It distinguishes the resulting operator algebra from the source's componentwise product, and states the conditions needed for the relative metric and rotation flow.

## R4 source connection

[R4](../03-lambda-reference/NATIVE_BOND_BRIDGE_R4.md) gives the exact coordinate map from the existing native EMK pair to R3's seam coordinates, including the relative metric. It distinguishes the represented quarter-turn from the central cut-complex scalar in the full native matrix algebra.

## R5 observable depth coordinates

[R5](../03-lambda-reference/PROFILE_RECOVERY_R5.md) identifies the finite list of positive depth products encoded by a calibrated boundary response jet. Its coordinate data are the odd coefficients at zero coupling. Exactly one additional coefficient resolves one additional cell; no change to the native state algebra or the source manuscript's unresolved field construction is claimed.

## R6 uncertain depth coordinates

[R6](../03-lambda-reference/INTERVAL_PROFILE_RECOVERY_R6.md) replaces exact coefficient coordinates by rational enclosures and records which positive product coordinates remain identifiable at the supplied precision. No additional primitive state algebra is introduced.
