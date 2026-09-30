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

Choose a concrete observable pair and derive its response from the Way-4 generator. Then determine whether the fixed reference of Way-3 constrains the response through an invariant dimensionless relation.

Role in the four-route program: connect an algebraic or operator result to a quantity that can be compared across states or measured. See [the shared assessment](../ASSESSMENT.md).

## R1 development

[RESPONSE_RETURNS_R1.md](RESPONSE_RETURNS_R1.md) identifies closed-return products and equal-endpoint path ratios that are independent of local references. It also explains why ordinary derivatives along one common tangent give identity returns and why conservation is an additional condition.
