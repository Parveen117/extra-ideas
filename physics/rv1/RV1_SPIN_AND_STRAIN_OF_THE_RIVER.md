# RV1 — the river of the frame has a strain and a spin: the law reads the strain, loops and gyroscopes read the spin

Stage of the physics line. Symbolic algebra (sympy): `rv1_spin_and_strain_of_the_river.py`, `test_rv1.py` (5 tests).

The owner's statement: Coriolis is a coupling of the radial and the tangential; a new dimension coupled
with the frame of measurement produces curvature — native or local.

## Sources used (read, unchanged)

| Source | Statement used |
|---|---|
| this line, TP1 | law Q = ¼I₁ + ½I₂ − I₃ on the frame's order defect |
| MO1, CO1 | frame class: flat slices, one time, e₀ = ∂_t + v·∂_X; count dτ² = dt² − \|dX − v dt\|²; lapse equation |
| CV1 | connection of the frame from its order defect |
| RMG9 (response-geometry series) T1, T3 | L = m(1 + uK + vS) + wR; the quadratic form is blind to w; a closed protocol returns 2w·Area |
| DM1, GR2 | the third cut; the gyroscope |

## Statements

Write ∂_i v_j = S_ij + W_ij: strain (symmetric: expansion and shear) and spin (antisymmetric: ½ curl v).

**T1 (the law reads the strain only).** [e₀, e_i] = −(∂_i v_j) e_j, and for every river v(t, X)

  Q = S:S − (tr S)².

The spin does not appear. (I₁, I₂, I₃ each contain it; only the combination selected by equivalence in TP1 drops it.)

**T2 (a rigid turn is silent in the law).** v = Ω ẑ×X: I₁ = 4Ω², I₂ = −2Ω², I₃ = 0, Q = 0.
Adding a rigid turn to any river leaves Q unchanged.

**T3 (what a free reading feels).** For a slow reading of largest count (MO1),

  acceleration = ∂_t v + ∇(½v²) + (curl v) × velocity.

One formula holds both pulls: for the fall v² = r_s/r it gives r_s/2r² inward; for the rigid turn it gives
Ω²ρ outward. The sideways push is the Coriolis term, with field curl v = 2Ω.

**T4 (gyroscope).** The frame's own connection turns the triad along e₀ at exactly the spin W.

**T5 (loop).** Around a circle in a rigidly turning river the two-way time difference is exactly
2·flux(curl v)/(1 − β²): RMG9-T3's 2w·Area with w = Ω, divided by the squared clock factor.

**T6 (the law's equation for the river).** Varying v inside the class gives curl curl v = 0.
Every irrotational river satisfies it. For an axial turn ω(r) it leaves two cases: ω constant (the rigid
turn) and ω ∝ r⁻³.

**T7 (the r⁻³ river).** Its Coriolis field is a dipole, 2J(3z X − r² ẑ)/r⁵, with no divergence and no
curl. Unlike the rigid turn it has strain: S:S = 18 J² sin²θ / r⁶.

**T8 (refusal at second order).** With radial fall b(r) and this river together,
Q = −(2/r²)(r b²)′ + 18 J² sin²θ/r⁶. No b(r) makes it vanish at every angle. Flat slices with one time —
exact for the point mass (MO1) and for the expanding space (CO1) — cannot hold a turning source beyond first order.

**T9 (the plane is RMG9's element).** For v = L X with L = m(1 + uK + vS) + wR:
expansion = 2m, curl = 2w, and Q = −2 det(symmetric part) = −2m²(1 − u² − v²). The R generator is the
Coriolis coupling of the two directions of the plane; its axis is the direction out of the plane (the
third cut of DM1).

## Answer to "native or local"

Both, and they separate exactly. The spin is local: it is a curvature of the river (curl v ≠ 0), every loop
and every gyroscope reads it, and the law does not. A Coriolis coupling by itself — a rigid turn — makes
no gravity. The strain is native: it is all the law reads. The radial–tangential coupling becomes native
curvature when the turn rate changes with radius, and then the law itself fixes the change as r⁻³.

## Numbers (illustration)

```text
gyroscope in a polar orbit at 7020 km, Earth's turn (J = 5.86e33 kg m²/s)   40.9 mas per year
two-way lag around 1 m² at the pole                                          3.2e-21 s
```

The first is the frame-dragging part only (GR2 gave the other part); the published expectation for the
orbiting-gyroscope experiment was 39 mas per year (general knowledge; orbit and J differ slightly).

## What is put in, what is not claimed

* This is the known physics of a turning frame and a turning source (Coriolis and centrifugal terms,
  loop lag, dragging of frames). Nothing new is predicted.
* T6 is the variation of the law with respect to the river inside the frame class — one of the law's
  equations, not all of them. That the constant in ω = 2J/r³ is the source's angular momentum (with G/c²)
  is taken from general knowledge, not derived.
* T8 is shown for b(r) with an axial ω(r) only. What replaces flat slices at second order is not built.
* "The axis is a new dimension" is the plain fact that the curl of a planar river points out of the plane.

## Open gate

Second order: let the triad strain as well (a unit block on the slices, GB1) and test whether the law then
closes for a turning source.
