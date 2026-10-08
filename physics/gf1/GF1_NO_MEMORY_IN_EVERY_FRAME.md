# GF1 — The two no-memory conditions hold for every frame: one plane leaves no memory, a turn rate leaves no memory, and nothing else is needed

Monty Dabas. 8 October 2026. Python 3.12. Symbolic algebra (sympy), with the functions of CV1 and TP1.

MM1 obtained TP1's coefficients 1 : 2 : −4 from two no-memory conditions, on flat slices with one time. Its claim
boundary lists "the same for general frames" as not shown. This stage shows it.

Sources read before building: **TP1** (the three invariants of the order defect; T2: equivalence selects
1 : 2 : −4; Q = −R + boundary), **MM1** (U6: the two conditions on flat slices), **SE1**, **MA1** (the pole boundary
term of the spherical frame), **PR3** (the accelerated frame), **SW1-T2** (constant swirl), **GE2-T5**, **IN1-T2**,
**EMK-1** (determinant in a plane), **FR1** (three cuts).

## The conditions are algebraic

Q = a₁I₁ + a₂I₂ + a₃I₃ is a quadratic form on the order defect c_ab^k at a place — 24 numbers, for any frame at
all. So the two conditions can be put to all frames at once.

**P — one plane.** Let the order defect lie in one plane: only [e_a, e_b] = c_ab^a e_a + c_ab^b e_b, for one pair
(a, b). For each of the six planes

```text
Q = ( 2a₁ + a₂ + a₃ ) × ( a non-zero form ) .
```

This covers a frame accelerated along one cut, a stretch along one cut, and a frame turning inside one plane of
space from place to place.

**T — a turn rate.** Let the cuts turn along e₀ at any rates: c_0i^j = −c_0j^i. Then Q = (2a₁ − a₂)·2(w₁² + w₂² + w₃²).

**PT.** "An order defect in one plane leaves no memory" and "a turn rate leaves no memory" hold together exactly
for a₁ : a₂ : a₃ = 1 : 2 : −4. No slices and no special frame are assumed.

## What the selected law is made of

**S.** In TP1's law no in-plane part of the order defect appears squared, and two in-plane parts multiply each
other only when their planes share exactly one direction. The memory of a frame is always *between* two planes;
one plane alone has none. On flat slices this is SE1's −2·Σ_{i<j} K_ii K_jj.

## Frames of flat space

```text
uniformly accelerated frame (time, x)        one plane        Q = 0 at every place          F
frame turning in time                        turn rate        Q = 0 at every place          F
polar frame in one plane of space            one plane        Q = 0 at every place          F
spherical frame                              two planes       Q ≠ 0, a pure boundary term:  Q = −2·div(trace vector)
```

For each of the first three a law with other coefficients gives a non-zero density. The last line is MA1's pole
term and TP1-T2's "boundary term": when the frame change spreads over two planes the memory density is not zero
place by place, but it is a flux and adds to nothing.

## What this settles

```text
MM1: "general frames — not shown"         the two conditions are algebraic and select 1 : 2 : −4 for every frame
equivalence of local frames (TP1-T2)      and the two no-memory conditions give the same law; the premises of the
                                          line reduce to one kind: no memory in a shared record
```

## What is put in

- TP1's starting point: the law is coordinate-free and quadratic in the order defect (three invariants).
- The reading of the two cases as records without memory: a turn rate as a deterministic record (GE2-T5), an order
  defect in one plane as a single reading (IN1-T2). The algebra does the rest.

## What is not shown

- That gravity has no content in two dimensions, and that this law is the only one with that property among the
  three-term laws, is known (general knowledge). The line's part is that the selection is stated as no-memory
  conditions of its own record rule, and S.
- Why the law is quadratic in the order defect is not derived.
- For frame changes that spread over more than one plane the statement is TP1-T2's (boundary term), not a
  place-by-place zero.

## Claim boundary

```text
ONE-PLANE ORDER DEFECT: Q ∝ 2a₁ + a₂ + a₃ , ALL SIX PLANES                                PROVED (symbolic)
TURN RATE: Q ∝ 2a₁ − a₂                                                                    PROVED
THE TWO CONDITIONS ⇔ 1 : 2 : −4 , FOR EVERY FRAME                                          PROVED
NO SQUARES OF IN-PLANE PARTS ; PRODUCTS ONLY BETWEEN PLANES SHARING ONE DIRECTION           PROVED
ACCELERATED, TURNING AND ONE-PLANE POLAR FRAMES: Q = 0 ; SPHERICAL FRAME: BOUNDARY TERM     PROVED (symbolic)
WHY QUADRATIC                                                                              NOT DERIVED
THE RECORD RULE                                                                            PREMISE
```

## Reproduce

```text
python gf1_no_memory_in_every_frame.py
python -m unittest test_gf1
```
