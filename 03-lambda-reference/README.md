# Way-3: A constant reference lambda

Source: Reference-Scale Normalization in [4ways.tex](../4ways.tex). [source-excerpt.tex](source-excerpt.tex) preserves that section verbatim as a LaTeX fragment.

## Core idea

Hold a reference `lambda_*` fixed and express states or quantities through coefficients and scale orders:

`X_i = n_i lambda_*^{k_i}`.

This is the route closest to the intuition that lambda is a constant underneath changing descriptions.

## Research judgment

This is an ambitious foundational direction. Its decisive question is whether the framework selects an invariant lambda, rather than allowing the reference and coefficients to be adjusted together. The supplied manuscript declares the reference but does not select its value.

The fixed reference and a changing response `chi_XY` can coexist. A relation between them must be derived; the shared use of the letter lambda does not identify them.

## Two distinct mathematical models

| Model | Necessary structure | Remaining freedom |
| --- | --- | --- |
| Scalar normalization | Nonzero reference and admissible coefficients | Free coefficients may absorb any allowed reference rescaling |
| Discrete valuation | A valued field, valuation-zero units, and a chosen uniformizer | Multiplication of the uniformizer by a unit preserves the scale-order structure |

Do not identify these models without a map and its hypotheses. If coefficients are restricted to integers, a lattice, or a compact unit group, state those restrictions explicitly and determine the surviving reference freedom.

## Next target

Use [IDENTIFIABILITY.md](IDENTIFIABILITY.md) as the first proof obligation. Derive a native relation that selects a dimensionless invariant, establish its uniqueness, and connect it to a response observable. A spectral ratio or return invariant is a candidate class of objects to investigate, not a value already obtained by this draft.

Role in the four-route program: the main question about scale selection. Way-4 supplies tools for addressing it and Way-2 supplies its response interpretation.
