# RC1 — Rational is closure, irrational is the turn left over

Monty Dabas. 8 October 2026. Python 3.12. Exact rational arithmetic; sympy for identities with symbols.

The owner's statement: in the rational there is exact closure, commensurability; in the irrational, a turn that
remains, non-closure. A point has no dimension; a circle is curved in one dimension and needs one more to be
laid flat; what is drawn as a circle on a flat surface means that surface has an imaginary part; in the same way
the ellipse will be the shortest path one dimension higher, in four.

Sources read before building: **MO1-M4/M5** (circles; near a circle (in-out rate / round rate)² = 1 − 3·r_s/r;
last stable circle 3·r_s; light circle (3/2)·r_s), **OA1-T3** (carried direction: count factor² = 1 − (3/2)·r_s/r),
**PT1** (exact turns of the rational block: quarter turns × one turn per prime p = x² + y²; **P4: no exact turn
repeats**), **QC1-T3** (a content is silent when q·Θ is a whole number of turns), **SW2** (the centre displaced by
ι·a; R = r − ι·a·cos θ; the ring), **DG1** (the diagonal), **PR1** (response space = velocity space).

## 1. Closure on the line's own circles

**R1 — two ladders.** With x = r_s/r:

```text
orbit              ( in-out rate / round rate )²            = 1 − 3x          MO1-M5
carried direction  ( turn against the round basis / 2π )²   = 1 − (3/2)·x      OA1-T3
```

The orbit closes exactly when the first ratio is a rational p/q; then r/r_s = 3q²/(q² − p²) and it closes after q
circuits. At any other radius it never closes: per in-out period the orbit goes round by 2π·(round rate)/(in-out
rate), an irrational number of turns, and the part beyond one turn is MO1's advance. The carried direction returns
exactly when the second ratio is rational: r/r_s = 3q²/(2(q² − p²)).

**R2 — the zero ends.** p = 0 on the first ladder is r = 3·r_s, MO1's last stable circle. p = 0 on the second is
r = (3/2)·r_s, MO1's light circle. The two special circles of the field are the zero points of the two ladders.

**R3 — both at once.** Orbit and carried direction close on the same circle exactly when

```text
1 − 3x = (P/Q)²   and   1 − (3/2)x = (C/Q)²      ⇔      P² + Q² = 2C²
                                                 ⇔      A² + B² = C² ,   A = (Q+P)/2 ,  B = (Q−P)/2 .
```

These are the right triangles with whole sides, that is, the exact turns u = (A + B·R)/C of PT1. So the doubly
closed circles are indexed by PT1's exact turns, one family per prime p ≡ 1 mod 4:

```text
u₅  = (3 + 4R)/5       r = (49/16)·r_s        orbit 1/7 ,  direction 5/7
u₁₃ = (5 + 12R)/13     r = (289/80)·r_s       orbit 7/17 , direction 13/17
u₁₇ = (15 + 8R)/17     r = (529/160)·r_s      orbit 7/23 , direction 17/23
```

and the orbit's ratio is tan(45° − φ), φ the angle of the exact turn: the departure from the diagonal of DG1.

**PT1-P4, recalled.** A rational point of the circle that is not a quarter turn never returns. "Rational" has two
meanings and they pull opposite ways: a rational *ratio of rates* closes; a rational *point* is an irrational
turn. R3 is where the two meet.

## 2. A point, a circle, an ellipse

**R5 — a point displaced along ι is read as a sphere.** In n dimensions the distance from a point moved by ι·a
along one axis vanishes on a sphere S^(n−2) of radius a in the plane across that axis: two points for n = 2, a
circle for n = 3 (SW2's ring), a sphere for n = 4. A circle on real flat space is a point whose place has an
imaginary part.

**R6 — the ellipse.** The curves on which the real part of that distance is constant are the ellipses
X²/(r² + a²) + Z²/r² = 1, and their foci are the two read points of R5. The circle is the case a = 0. An ellipse is
a circle about a point that has an imaginary part.

**R7 — the ellipse as a shortest path in four dimensions.** At first order MO1's orbits are ellipses (M5 without
its last term). Their velocities, lifted from the three-dimensional velocity space to the 3-sphere in four
dimensions, run along great circles: the shortest paths there. Every great circle closes: this is the ratio 1 of
R1 at x → 0. The last term of M5 is what makes the ratio √(1 − 3x) and leaves a turn over.

## 3. What this gives

```text
closure            a rational ratio of two rates of one history                       R1
left-over turn     the irrational part; MO1's advance, OA1's turn                     R1
special circles    the zero points of the ladders                                     R2
exact turns        index the circles where everything closes                          R3
circle, ellipse    a point with an imaginary part, read on real flat space            R5, R6
ellipse            shortest path one dimension higher                                 R7
```

## What is put in

- MO1's and OA1's formulas (for a centre at rest), PT1's classification, SW2's displaced distance.

## What is not shown

- That nature prefers the rational radii. That is tested against measured rates in QP1, not derived.
- "A point with no dimension is pure irrational; observed in one dimension it is rational" is not formalised here.
- R5–R7 are elementary or known geometry (the lift of R7 is known; general knowledge, no source re-read). What is
  the line's own is R2 and R3.
- The doubly closed circles are not known to be observable.

## Claim boundary

```text
ORBIT CLOSES ⇔ √(1 − 3x) RATIONAL ; DIRECTION RETURNS ⇔ √(1 − 3x/2) RATIONAL            PROVED (exact)
ZERO ENDS = LAST STABLE CIRCLE AND LIGHT CIRCLE                                         PROVED (exact)
DOUBLY CLOSED CIRCLES ⇔ EXACT TURNS OF PT1 ; ORBIT RATIO tan(45° − φ)                   PROVED (exact, both directions, checked to Q < 200)
ι-DISPLACED POINT READ AS S^(n−2) ; ELLIPSE = CONSTANT REAL PART OF ITS DISTANCE        PROVED (symbolic)
FIRST-ORDER ORBITS = GREAT CIRCLES OF THE 3-SPHERE                                      PROVED (symbolic; known)
NATURE SELECTS THE RATIONAL POINTS                                                      NOT DERIVED — see QP1
```

## Reproduce

```text
python rc1_rational_closure.py
python -m unittest test_rc1
```
