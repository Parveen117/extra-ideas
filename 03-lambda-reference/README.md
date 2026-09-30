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

## R1 development

[RETURN_INVARIANTS_R1.md](RETURN_INVARIANTS_R1.md) now constructs reference-independent quantities `Lambda_C` attached to specified return processes, classifies the finite scalar model, and proves a conservation criterion. These are computable return multipliers; a law selecting a universal `lambda_*` is still required.

## R2 development

[AGHORA_RETURN_R2.md](AGHORA_RETURN_R2.md) constructs a specified return family from the source Aghora relations. Its eigenvalues satisfy `r^2=1`, while a bond projection can produce a variable effective coefficient. The remaining selection question is to derive the native retained subspace and its relationship to the return, rather than choose one sign or projection by convention.

## R3 development

[SEAM_BOND_COMPLEX_STRUCTURE_R3.md](SEAM_BOND_COMPLEX_STRUCTURE_R3.md) selects a bond under the added seam-exchange rule, proves `[A,B]^2=-I`, and derives a relative metric and the generator form `G=omega[A,B]` under separate compatibility conditions. The rate, phase, sign choice, and universal reference remain unselected. The next question is whether the original operations force these added conditions.

## R4 source connection

[NATIVE_BOND_BRIDGE_R4.md](NATIVE_BOND_BRIDGE_R4.md) realizes the bond rule using the existing native exact-cut response. The removed space follows within that source family. Source selection, the coupling profile, and dynamical calibration remain open; the earlier iota and cut-synthesis results are explicitly credited.

## R5 identification and source-selection limit

[PROFILE_RECOVERY_R5.md](PROFILE_RECOVERY_R5.md) gives a constructive inverse problem for arbitrary positive depth dependence. Finite calibrated response data determine a precise prefix. It also constructs, for every finite prefix length, distinct bounded positive sources with identical available coefficients and the same exact closed bond. A selected bond and a successful finite fit therefore do not select the full source in this observation protocol. The additional probe and tail error budget make that distinction testable within the native model; a universal `lambda_*` remains unselected.

## R6 a certified recovery boundary

[INTERVAL_PROFILE_RECOVERY_R6.md](INTERVAL_PROFILE_RECOVERY_R6.md) determines how far uncertain response data certify a positive depth profile. It distinguishes an unresolved next product from an incompatible source model, rather than treating a fitted central value as selection. The universal reference still needs an independent source law.
