# QC1 — Integers come from the reading: content, closure and the weight ladder

Monty Dabas. 7 October 2026. Python 3.12. First stage of the
quantum–classical line; continues PH3.

**Question put to this stage.** If a shared evolution must be recordable in
a pure reading, the return has to be invisible to it: Θ ∈ 2πℤ in some
sense. Do the integers that appear match the discrete series of SL(2,ℝ),
the group of the response space (RMG10)?

**Answer.** Yes for the ladder, with a correction to where the integers
come from. The response plane is contractible, so the geometry imposes no
integer at all. Every integer below enters through the **content q of the
reading** — the degree at which the probe is read.

## 1. The two readings are a canonical pair

Offsets w = (δs, δv) and p = H·w = (δT, −δP).

**T1.** The E reading moves p by dH·w and fixes w; the M reading moves w by
d(H⁻¹)·p and fixes p. Their equal share keeps the pair on the equilibrium
graph p = H(t)w, acts there as the native transport dw = −½H⁻¹dH·w, and
conserves the pairing p·w.

**T2 (the return is their bracket).** With A_i = ∂_iH and B_i = ∂_i(H⁻¹),

```text
F_½ = ¼ H⁻¹ (A_s B_v − A_v B_s) H .
```

The curvature that produces Θ is the bracket of the generators of the two
readings, ½wᵀA w and ½pᵀB p. Two readings that commute give no return.

## 2. Content

**T3.** Let z = y₁ + iy₂ in the H-orthonormal probe plane. Under a return
of angle Θ the monomial z^a z̄^b is multiplied by e^{i(a−b)Θ}. A reading of
degree m therefore splits into contents q = a − b ∈ {−m, −m+2, …, m}, and

```text
content q is silent  ⇔  qΘ ∈ 2πℤ .
```

- q = 0 is always silent: the second-order energy ½wᵀHw (PH3).
- First moments (q = ±1) are silent only at Θ ∈ 2πℤ.
- The covariance of the probe (q = 0, ±2) is silent at Θ ∈ πℤ.

## 3. The ladder

On a circle about isotropy at hyperbolic radius d (RMG3),
Θ = −π(cosh d − 1). Put k = q/2.

**T4.** Content q is silent exactly on the circles

```text
cosh d_n = 1 + 2n/q ,        k·cosh d_n = k + n ,        n = 1, 2, 3, …
```

So on the hyperboloid of radius k the silent circles sit at the heights
k + n. These are the weights of the discrete series of SL(2,ℝ) with index
k (the comparison target; named classical input, not used in any proof).
Odd contents give half-integer k: the double cover.

```text
q = 1, k = ½ :  heights 3/2, 5/2, 7/2, …        q = 2, k = 1 :  heights 2, 3, 4, …
q = 3, k = 3/2: heights 5/2, 7/2, 9/2, …        q = 4, k = 2 :  heights 3, 4, 5, …
```

For the covariance (q = 2) the first silent circle has Klein radius²
3/4, an eigenvalue ratio of 7 + 4√3 ≈ 13.9.

## 4. A real fluid: the first covariance-silent CO₂ cycle

The PH1 ellipse with its lowest temperature at T_c + δ:

```text
δ [K]     5.87     3       1       0.5     0.25    0.1     0.05    0.03    0.02    0.01    0.001
Θ        −0.926  −1.261  −1.782  −2.083  −2.363  −2.697  −2.917  −3.057  −3.152  −3.281  −3.448
```

Θ crosses −π at **δ = 0.0210 K**. On that cycle the covariance of the probe
returns unchanged (relative change 3·10⁻⁶) while every first moment
returns reversed (w → −w).

The reference equation of state is not reliable to 0.02 K at the critical
point; the number shows that a silent cycle exists in this family and
roughly where, not a measured value.

## 5. Certificate

Exact rational arithmetic: T1 and T2 at three points of a quartic convex
energy (graph preservation, pairing, bracket = curvature, non-commuting
witness); T3 as polynomial identities over ℚ(i) for degrees 1–5 with
rational rotations; T4 as rational identities for q = 1…4. Five tests,
two tampers (a wrong share leaves the graph; a non-rotation breaks the
content law).

## 6. Claim boundary

```text
THE TWO READINGS ARE A CANONICAL PAIR; RETURN = THEIR BRACKET          PROVED
CONTENT LAW; SILENCE ⇔ qΘ ∈ 2πℤ                                        PROVED
SILENT CIRCLES AT HEIGHTS k + n, k = q/2                               PROVED; matches the discrete-series weights
INTEGERS FORCED BY THE GEOMETRY ALONE                                  NO — the plane is contractible; they come from q
A UNIT OF RETURN (ħ, k_B) DERIVED                                      NO — nothing here has a scale
COVARIANCE-SILENT CO₂ CYCLE                                            COMPUTED (δ ≈ 0.021 K, equation of state not reliable there)
```

What this stage supports: quantization as a property of the reading, not
of the geometry. What it does not reach: why a physical reading has
integer content with a fixed unit. That needs the pair (w, p) to carry a
scale — the fluctuation relation ⟨δx_i δy_j⟩ = k_BT δ_ij on the thermal
side — and is the next gate.

## 7. Reproduce

```text
python -m unittest test_qc1_exact
pip install CoolProp numpy && python qc1_reading_content_closure.py
```
