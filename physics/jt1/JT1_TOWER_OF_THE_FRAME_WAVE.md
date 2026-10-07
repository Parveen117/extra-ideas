# JT1 — The tower of the frame wave: flat where the axis is rigid, curved where it turns

Monty Dabas. 8 October 2026. Python 3.12. Symbolic algebra (sympy);
exact rationals for the tower observer.

GW1 stopped at second order and said what lies beyond was not built.
The owner's reply: the higher orders are the Jacobian tower, and that is
where the framework's curvature and flatness decide.

Sources read before building: RKF **Theorem 42** (cut-graded
λ-Jacobian tower: L_{n+1} = ∇L_n; termination for polynomial data,
eq. 3.4; generating object Σ t^k/k! ∇^k, eq. 3.5; tower observers and
target repair, Theorem 7.1; higher-layer recovery, Prop. 8.1);
**NT-3** (X = G⁻¹δG); **RMG1** (response element ↔ point of the
hyperbolic plane; curvature needs the principal axis to turn; link to
Hessian-power rigidity, [H, ∂H] = 0); **TP1**, **GW1**, **GR2-G1**.

## Results

Transverse frame block M(t, z); response element G = MᵀM.

**J1 — the law, exactly, at every order.**

```text
Q = ¼ [ tr(X_t²) − tr(X_z²) ]  −  [ (∂_t ln det M)² − (∂_z ln det M)² ] ,      X = G⁻¹ ∂G .
```

The law of TP1 on a travelling frame is the invariant square of NT-3's
X — the same X whose commutator is the native curvature. Verified
against the frame computation on general frames, volume term included.

**J2 — closed form.** For unit volume, G = Exp(−σ(cos φ K + sin φ L)):

```text
Q = ½ [ σ_t² − σ_z² + sinh²σ ( φ_t² − φ_z² ) ] .
```

σ is the size of the wave (the anisotropy ℓ of RMG1), φ the direction
of its principal axis.

**J3 — flat tower.** If the axis does not turn — [G, ∂G] = 0, the
rigidity condition — the law is exactly ½(σ_t² − σ_z²). Every layer of
the tower above the second vanishes (Theorem 42, eq. 3.4, in the Exp
chart). A wave of one polarisation never acts on itself.

**J4 — curved tower.** If the axis turns,

```text
sinh²σ = σ² + ⅓ σ⁴ + (2/45) σ⁶ + … ,
```

so the first layer beyond GW1 is ⅙ σ⁴(φ_t² − φ_z²), and every higher
layer comes from the one closed form (the generating object of
eq. 3.5). The interaction is the curvature of the response plane and
nothing else.

**J5 — one-way waves.** For any σ, φ depending on t − z only, Q = 0 at
every order.

**J6 — the tower observer and equivalence.** At a point of the static
field, with what a change of frame can add,
m = r_s/r + c₀ + c₁(r − r₀):

```text
layer 1 (value)             blind to r_s
layer 2 (first derivative)  blind to r_s
layer 3 (second derivative) recognises it:   r_s = m″ r₀³ / 2
```

GR2-G1 (memory removed at a point) and CV1 (curvature not removed) are
Theorem 42's target repair: the source is first recognised at the third
layer.

## Reading

"Nonlinear" has an exact meaning here: the turning of the principal
axis of the response element. Flatness of the tower is rigidity;
its curvature is the curvature of the response plane of RMG1. The
question left open in GW1 — what happens when counts meet — has this
answer at the level of the law: waves of one polarisation, and waves
travelling one way, do not interact at any order; interaction needs
two polarisations meeting, and its strength is fixed by sinh²σ.

## What is not shown

- The frame family has only the transverse block varying; the other
  components of the frame, which carry the conditions the full law
  imposes on such waves, are held fixed. J1–J5 are statements about the
  law on this family, not complete solutions.
- These features — no self-action of one polarisation or of one-way
  waves, and a hyperbolic-plane form for two — are known for plane
  gravity waves. The framework's own part is the identification with
  NT-3's X, RMG1's plane and Theorem 42's tower.
- What interacting counts do (the quantum side of J4) is not built.

## Claim boundary

```text
Q = ¼[tr X_t² − tr X_z²] − VOLUME TERM, ALL ORDERS                           PROVED (symbolic, general frames of the family)
CLOSED FORM WITH sinh²σ                                                      PROVED
RIGID AXIS ⇔ TOWER STOPS AT SECOND ORDER                                     PROVED
FIRST HIGHER LAYER ⅙σ⁴(φ_t² − φ_z²); ALL LAYERS FROM ONE FORM                 PROVED
ONE-WAY WAVES: Q = 0 AT EVERY ORDER                                          PROVED
SOURCE RECOGNISED FIRST AT LAYER THREE                                       PROVED (exact)
FULL FRAME (ALL COMPONENTS); INTERACTING COUNTS                              NOT BUILT
A PREDICTION BEYOND THE KNOWN LAW                                            NONE
```

## Reproduce

```text
pip install sympy
python jt1_tower_of_the_frame_wave.py
python -m unittest test_jt1
```
