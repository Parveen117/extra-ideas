# CD1 — Closure picks three: the least-cost memory closes a bound history only in three dimensions

Monty Dabas. 8 October 2026. Python 3.12. Symbolic algebra (sympy).

The line carries "three space dimensions" as an input (MC1's claim boundary) and as an identification with the
number of cuts (FR1's). RC1 states the owner's closure rule: rational ratio of rates = closure. This stage asks
what that rule says about the number of dimensions.

Sources read before building: **MC1-T1** (least-cost reading in d dimensions: r^−(d−2); log r for d = 2; claim
boundary: three dimensions is an input), **FR1** (three cuts and no fourth; "dimension = number of cuts" is an
identification, not derived), **MO1-M1/M3/M5** (static form, radial law, orbit equation), **GR1** (N² = 1 − m),
**RC1** (closure), **QC1-T3** (silence in content 1: a return of a whole number of turns), **OA1** (count factor
on a free circle), **CV1** (the frame law).

## Results

**D1, D2 — any memory.** With MO1's form and any memory m(u), u = 1/r, a free circle has

```text
( in-out rate / round rate )²  =  1 − m − 2u·m′ − u·(1 − m)·m″ / m′ .
```

**D3 — the least-cost memory of MC1.** For m = μ·u^(d−2):

```text
( in-out rate / round rate )²  =  ( 4 − d ) − d·m            (d ≥ 3) ,
                                  2 − 2m − 2μ                (d = 2, m = c + μ log u) .
```

In d = 3 this is MO1-M5's 1 − 3m.

**D4 — closure at first order.** Far from the centre (m → 0) the ratio is √(4 − d):

```text
d = 2      √2          irrational: no bound history near a circle ever closes
d = 3      1           every bound history closes, after one circuit
d = 4      0           no in-out motion; with any memory the ratio² is −4m: no stable circle at any radius
d ≥ 5      none        ratio² negative: no stable circle
```

So among all dimensions, the least-cost memory gives closed bound histories only in three. In QC1's words: a
bound history is silent in content 1 (returns after one circuit) only in three dimensions.

**D5 — the special shells.** On a free circle the count factor² is 1 − (d/2)·m (OA1-T3 in d dimensions). The last
stable circle, the light circle and the shell where a held clock stops are at

```text
m = (4 − d)/d ,      2/d ,      d/d :        equal steps of (d − 2)/d in the memory, in every dimension.
```

In three dimensions they are 1/3, 2/3, 1. The first of them exists (positive memory) only for d = 3.

**D6 — the same exponent from the frame law.** For d = 3, 4, 5 the static form with memory B·r^−(d−2) has zero
contracted curvature, and with another exponent it does not (coordinate route). Least cost (MC1) and the frame
law (CV1) give the same memory in each dimension checked.

## What this settles

```text
MC1:  "three space dimensions"            INPUT        →   the one dimension in which the least-cost memory closes a bound history
FR1:  "dimension = number of cuts"        IDENTIFICATION →  a second, independent reason for three: closure
```

It is conditional on the closure rule (RC1) being a requirement on records, which QP1 tests and nothing yet
derives.

## What is put in

- MC1's least-cost memory and its assumptions; MO1's form used with a general memory; the closure rule.

## What is not shown

- That bound orbits close only for the 1/r law, and that stable orbits need at most three dimensions, are known
  results (general knowledge; no source re-read). The line's part is that its own least-cost memory and its own
  closure rule meet there, and D5.
- First order only for the closure statement: with the memory in, three dimensions closes on the rungs of RC1.
- D6 is checked for d = 3, 4, 5 by the coordinate route; CV1's frame route is written for three dimensions.
- Nothing here says why records must close.

## Claim boundary

```text
RATIO² FOR ANY MEMORY ; (4 − d) − d·m FOR THE LEAST-COST MEMORY                      PROVED (symbolic, any d)
FIRST-ORDER CLOSURE ONLY IN d = 3 ; √2 IN d = 2 ; NO STABLE CIRCLE FOR d ≥ 4          PROVED
SPECIAL SHELLS AT (4−d)/d, 2/d, 1 : EQUAL STEPS ; THIRDS IN THREE DIMENSIONS          PROVED
LEAST COST AND THE FRAME LAW GIVE THE SAME EXPONENT                                   CHECKED for d = 3, 4, 5
THREE DIMENSIONS DERIVED                                                              CONDITIONAL on the closure rule
THE CLOSURE RULE                                                                      NOT DERIVED
```

## Reproduce

```text
python cd1_closure_and_dimension.py
python -m unittest test_cd1
```

## Later note (SE1, 8 October)

The d = 2 row above uses MC1's least-cost memory, log r. SE1 shows that the frame law in two cuts leaves β
constant: no falling memory at all. Either way nothing closes in two dimensions; D6's agreement of least cost and
the frame law holds for d ≥ 3 only.

## Later note (FK1, 9 October)

D1 with the count factor² = 1 − d·m/2 gives (in-out/round)² − (seen − lost) = 3 − d: the law's F of a free circle
is its squared rate ratio in three cuts only. Nothing above is changed.
