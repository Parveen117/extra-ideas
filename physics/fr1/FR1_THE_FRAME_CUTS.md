# FR1 — The frame is what cuts: one invariant, many readings, and how many cuts a frame carries

Monty Dabas. 7 October 2026. Python 3.12, exact rational arithmetic only.
Continues MC1 and the synthesis.

Owner's statement for this stage: the geometry — circle, metric and the
rest — lives on the observer's frame, not on nature. The observer reads
one invariant through different cuts. If there is a frame, it is the
frame that cuts.

Sources read before building: **EMK-1** T1–T3, T7 (the block, its
determinant channels Δ∥, Δ⊥ defined against the seam reflection K; the
commutator with K as the seam-mixing detector); **EMK-T1** T1 (the
cut-swap is derived from the basis, not primitive); **SY1** S1, S5;
**PR1-T4** (propagation by a cut, mass by the turn); **PR2-T1/T4** (form
and sheets; the pair record); **RMG9-T3** (three rotation invariants);
the cut-complex unit ι (F00E: ι² = −1), the framework's primitive carrier.

## 1. One invariant, many readings

**T1.** Turn the cut inside the split plane, K_φ = cos φ·K + sin φ·RK.
The channels of a block read against K_φ are

```text
Δ∥(φ) = a² − b_φ² ,      Δ⊥(φ) = c² − d_φ² ,      b_φ = b cos φ + d sin φ ,   d_φ = −b sin φ + d cos φ ,
```

and their sum does not depend on φ. Δ∥ ranges over [a² − (b² + d²), a²].
What no turned cut can change: a, c and b² + d² (RMG9's three numbers).
Under every change of frame of a unit block, only E = a survives.

```text
block                    invariant     read along K      read across K
boost                    1             (1, 0)            (25/16, −9/16)
coin turn                1             (9/25, 16/25)     (9/25, 16/25)       the same against every turned cut
response                 23/9          (32/9, −1)        (3, −4/9)
```

The same boost is "all seam-compatible" for one cut and "partly
rotational" for another. A response's coupling r² = −Δ⊥/Δ∥ (RMG9) is
9/32 for one cut and 4/27 for another: it belongs to the cut. A turn is
read alike by all turned cuts, and differently by a boosted one (PR2).

**T2 (the geometry is the set of cuts).** Cuts are the traceless
elements of form +1 — one sheet; turns have form −1 — two sheets; state
tensors are null. The circle is the orbit of one cut under the turn; the
hyperbolic sheets are the orbits under boosts; the metric is the form
itself. None of it is a property of the block being read.

## 2. How many cuts a frame carries

A sector that propagates by cuts needs one cut per direction, all
anticommuting, so that the square is the sum of squares (PR1-T4).

**T3.**

```text
over the rationals:             two anticommuting cuts, K and RK; the only further anticommuting element is R, a turn
with the cut-complex unit ι:    three anticommuting cuts, K, RK and ιR; no fourth linear one
```

The block read with real coefficients has two directions and one turn.
Read on the framework's own cut-complex carrier, the turn becomes a third
cut. Three is the maximum.

**T4 (in three directions the mass is the other sheet).** With all three
cuts used for propagation, no linear element anticommutes with them: a
single block has no mass term. The conjugation J_c = (ιR)∘conj does —
it anticommutes with the three cuts, J_c² = −I, and

```text
( ξ₁K + ξ₂RK + ξ₃ιR + g J_c )² = ( ξ·ξ − g² ) I .
```

The elements anticommuting with all three cuts are exactly a
two-dimensional family, all of them conjugations. In three directions a
rest turn exists only as a coupling between the block and its conjugate
sheet — the pair of PR2.

## 3. Reading

- What is one is the scalar part. Speed against mass, seam against
  rotation, coupling, direction — these are how a cut divides it.
- Dimension is a property of the frame: the number of anticommuting cuts
  it carries. Two for a real reading, three for a cut-complex one.
- With three, MC1 gives the exponent d − 2 = 1: the 1/r of GR1 is the
  least-cost reading of a frame that carries three cuts.
- With three, there is no room left for a turn inside one block. Mass is
  then necessarily the meeting of a reading with its conjugate — the
  observer "sitting with the antimass".

## 4. Certificate

T1 on four blocks and five rational cut angles, with the components
computed two ways; the frame change by a boost. T2 on seven elements.
T3–T4 by exact nullspace computations: real 2×2 elements anticommuting
with K and RK; on the cut-complex carrier (as real 4×4) the linear,
all, and antilinear solutions (dimensions 0, 2, 2); the square law at
two rational points. Six tests; a non-cut and a non-anticommuting third
element are rejected.

## 5. Claim boundary

```text
CHANNELS DEPEND ON THE CUT; ONLY E SURVIVES EVERY FRAME CHANGE              PROVED
CUTS, TURNS, STATES AS FORM +1, −1, 0                                       PROVED
TWO CUTS (REAL), THREE (CUT-COMPLEX), NO FOURTH                             PROVED
NO LINEAR MASS WITH THREE CUTS; THE CONJUGATE MASS AND ITS SQUARE LAW       PROVED
"SPACE DIMENSION = NUMBER OF CUTS OF THE FRAME"                             IDENTIFICATION — fits PR1 and MC1; not derived
THAT PHYSICAL SPACE IS THE CUT-COMPLEX READING                              NOT CLAIMED
A SCALE; A VALUE OF ANY CONSTANT                                            NO
```

The algebra of three anticommuting cuts is a known one (general
knowledge). What is the framework's here is where it comes from: the
certified block, read on the cut-complex carrier, with the turn taken as
a cut.

## 6. Reproduce

```text
python fr1_the_frame_cuts.py
python -m unittest test_fr1_exact
```

## Later note (CD1, 8 October)

"Space dimension = number of cuts" is listed above as an identification. CD1 reaches three from another side:
the least-cost memory closes a bound history only in three dimensions. The two reasons are independent; that they
are one is not shown. Nothing above is changed.
