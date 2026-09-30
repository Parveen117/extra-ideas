# Way-1: Intrinsic multidimensional numbers

Source: the section of that name in [4ways.tex](../4ways.tex). [source-excerpt.tex](source-excerpt.tex) preserves it verbatim as a LaTeX fragment.

## Core idea

Represent a state by a vector in `F^n`, where `F` is an explicitly defined field. Use componentwise addition and multiplication. The coordinate idempotents `e_i` satisfy `e_i odot e_j=0` for `i != j` and decompose the algebra into its coordinate factors.

## What can be used now

Once the field is specified, this is a concrete commutative product algebra. It supports coordinate projections, matrix actions, and a precise starting state space. For `n>1` it has zero divisors; it is not itself a field.

## What needs development

The source's underlying field construction is incomplete, and its exact-shell grading fails. Begin with a defined discretely valued field and the filtration

`F^k(F^n) = {x : min_i v(x_i) >= k}`.

Construct the associated graded algebra through quotient spaces rather than declaring exact valuation shells to be graded pieces. Specify which changes of basis are admissible: arbitrary invertible linear maps need not preserve coordinatewise multiplication.

## Next target

Define the native state algebra and its admissible morphisms. Identify a structure beyond ordinary coordinate bookkeeping that is preserved by the maps to Way-2 and Way-4.

Role in the four-route program: a clear representation of states and components. See [the shared assessment](../ASSESSMENT.md) for explicit counterexamples and repairs.
