# MA1 — MO1's assumption follows from the line's own law and one condition far away

Monty Dabas. 8 October 2026. Python 3.12. Symbolic algebra (sympy), with the functions of CV1 and TP1.

The gravity thesis lists as its first open item: derive, or refuse, the assumption of MO1 — that the locally
pure frames share one flat space and one time. `UNUSED_RESULTS_MAP.md` proposed RKF theorum/31–32 for it.

Sources read before building: **MO1** §1 (the assumption), **GR2-G1** (the locally pure frame is boosted inward
with β² = m), **CV1** (curvature of a frame as order defect; contracted curvature zero in empty space),
**TP1** (the quadratic law selected by local frame equivalence, ¼I₁ + ½I₂ − I₃ = −R + boundary term, proved there
for the fall frame), **MC1** (m = r_s/r from least cost), **GB1**, **SC1**; RKF **theorum/31** (cut-variational
minimal observer) and **theorum/32** (cut-covariance event realisation).

On theorum/31–32: they fix the rank of a faithful observer and the covariance of a cut in the recognition
carrier. Neither speaks about how frames at different places fit together; they are not used here.

## The family

MO1's frame with the radial stretch left free:

```text
e⁰ = dt ,   e¹ = A(r)·( dr + b(r) dt ) ,   e² = r dθ ,   e³ = r sin θ dφ .
```

MO1 assumed A = 1 (flat slices). **T1:** this is no restriction on a static frame with clock factor N < 1 and
radial factor S: it is that frame with A = N·S and b² = (1 − N²)/(N·S)². So the assumption is the statement N·S = 1.

## Results

**T2 — TP1's identity holds on the whole family.** ¼I₁ + ½I₂ − I₃ = −R − 2·div(trace of the order defect),
pointwise, for every A and b. The selected law is the curvature law here too.

**T3 — the mixed component.** The contracted curvature has the mixed component −2·b·A′/(r·A²). In empty space
it vanishes, so wherever the frame falls (b ≠ 0), **A′ = 0**.

**T4 — the rest.** With A constant, every component vanishes exactly when A²·(r·b²)′ = 1 − A², i.e.
b² = (1 − A²)/A² + r_s/r.

**T5 — the condition far away.** If the frame is at rest far away (b → 0), then A = 1. The slices are flat, the
time is one, and b² = r_s/r.

**T6 — the other members.** A ≠ 1 is the same field read by frames that keep a speed far away,
b² → (1 − A²)/A². MO1's assumption selects, among them, the frames at rest far away.

## What this settles

```text
MO1's assumption (flat slices, one time)     =   N·S = 1                              T1
empty space, the line's law                   ⇒   A constant                           T2, T3
frames at rest far away                       ⇒   A = 1 and β² = r_s/r                 T4, T5
```

The assumption is not independent. Given the law of CV1/TP1, it follows for the static field of one centre from
the single condition that the pure frames are at rest far away; and the profile r_s/r, which MC1 obtains from
least cost, comes out of the same computation.

## What is not shown

- This is the known uniqueness of the static field of one centre, reached with the line's law. Nothing new is
  predicted.
- The law itself rests on TP1's premises (a quadratic, coordinate-free law; equivalence of local frames).
- Static and spherically symmetric only. Whether the pure frames of a general field fit together is not treated.
- Inside the horizon the falling form continues; the static form used in T1 does not (N < 1 is assumed there).

## Claim boundary

```text
ASSUMPTION OF MO1  ⇔  N·S = 1                                              PROVED
TP1's IDENTITY ON THE FAMILY WITH A FREE RADIAL STRETCH                    PROVED (symbolic)
EMPTY SPACE ⇒ A CONSTANT ; A²(r b²)′ = 1 − A²                              PROVED (symbolic)
AT REST FAR AWAY ⇒ A = 1 , β² = r_s/r                                      PROVED
THE ASSUMPTION FOR A GENERAL (NON-STATIC, NON-SPHERICAL) FIELD            NOT TREATED
RKF theorum/31–32 AS A ROUTE TO THIS                                       NOT APPLICABLE
```

## Reproduce

```text
python ma1_mo1_assumption_from_the_law.py
python -m unittest test_ma1
```
