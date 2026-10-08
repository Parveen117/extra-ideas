# CY1 — The cyclic form of three readings: one part is free of the reference, the other is not

Monty Dabas. 9 October 2026. Python 3.12, standard library only. Exact rational arithmetic.

Source object: the "cyclic tensor" of the catalog draft,

```text
T_AB = α_A β_B + β_A γ_B + γ_A α_B ,        α = ∇ln I ,  β = ∇ln J ,  γ = ∇ln K ,
```

its stated law "T_AB + T_BC + T_CA = 0 on the seam", and the tensor built from it,
D_AB = T_AB/(IJK) − ½ h_AB ln(IJK) + w² δ_AB.

Sources read before building: **R1** (return multipliers and ratios of equal-endpoint gains do not depend on
the reference), **Publications CID-1** (the split of a response depends on the metric, the obstruction does not).

## Results

**Y1 — relabelling.** A cyclic relabelling of (α, β, γ) leaves T unchanged; exchanging two of them gives the
transpose. The symmetric part is the same for all six orders; the alternating part changes sign with the order.

**Y2 — symmetric part.** With σ = α + β + γ,

```text
T + Tᵗ = σσᵗ − ααᵗ − ββᵗ − γγᵗ ,          2·tr_h T = h(σ,σ) − h(α,α) − h(β,β) − h(γ,γ) .
```

**Y3 — alternating part.**

```text
T − Tᵗ = (α − γ) ∧ (β − γ) .
```

It does not change when one covector is added to all three, and it is zero exactly when the two differences
are parallel.

**Y4 — readings.** For α = d ln I, β = d ln J, γ = d ln K:

```text
T − Tᵗ = d ln(I/K) ∧ d ln(J/K) .
```

A common factor of the three readings (I, J, K → cI, cJ, cK: a change of reference) leaves it unchanged. The
symmetric part T + Tᵗ changes by 2(σδᵗ + δσᵗ) + 6δδᵗ, δ = d ln c; this is zero only for δ = 0 or δ = −⅔σ,
that is, for the one factor c = (IJK)^(−2/3). For monomial readings in x, y the alternating part is the
determinant of the exponent differences times d ln x ∧ d ln y.

**Y5 — consequence for the draft.**

- The draft says T is not symmetric and gives its trace; both are right (Y2 contains the trace). It follows
  that D_AB as displayed is symmetric only where d ln(I/K) ∧ d ln(J/K) = 0, that is, where the three readings
  have at most one independent ratio. Elsewhere it cannot equal a symmetric curvature tensor.
- The displayed law "T_AB + T_BC + T_CA = 0 on the seam" is not an equation between objects of one type, and
  the draft does not say how the seam restricts it. Read on triples of index values it fails for general
  readings. What does vanish, and when, is Y3.

## Use

No stage uses this tool yet. What it offers: of the cyclic form of three readings, the part that no common
change of reference alters is one two-form, the area form of their two ratios (R1's invariants are likewise
ratios). The symmetric part depends on the reference, except for the one factor above. So T/(IJK) plus scalar
terms is in general neither symmetric nor reference-free. Compare CID-1, where the split depends on a choice
and the obstruction does not.

## What is not shown

- This is elementary multilinear algebra. Nothing here says which symmetric tensor, if any, should replace the
  draft's D_AB; that is a matter for the information–curvature source, which was not read for this stage.
- No statement about a seam, a metric h, or w is made beyond the trace formula.

## Claim boundary

```text
T + Tᵗ = σσᵗ − ααᵗ − ββᵗ − γγᵗ ;  T − Tᵗ = (α−γ)∧(β−γ)                                      PROVED
ALTERNATING PART = d ln(I/K) ∧ d ln(J/K), UNCHANGED BY A COMMON FACTOR                            PROVED
SYMMETRIC PART CHANGES BY 2(σδᵗ + δσᵗ) + 6δδᵗ ; ZERO ONLY FOR δ = 0 OR δ = −⅔σ                   PROVED
THE DRAFT'S CYCLIC LAW AS WRITTEN                                                                NOT KEPT
A REPLACEMENT FOR D_AB                                                                           NOT GIVEN
```

## Reproduce

```text
python cy1_three_reading_cyclic_form.py
python -m unittest test_cy1
```
