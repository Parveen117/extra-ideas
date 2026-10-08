# FK1 — The flip of a free circle is the squared ratio of its two rates; only in three cuts

Monty Dabas. 9 October 2026. Python 3.12. Symbolic algebra (sympy), numbers with mpmath.

ONE_LAW: S = R + D, F = R − D. For a free circle of a centre, seen R is the count factor² and lost D = 1 − R.
AC1-A3 found R = D at the last stable circle. This stage asks what F is on every circle, and what a turning
centre does to it.

Sources read before building: **T24-6.1** (S = R + D, F = R − D), **DO1**, **AC1-A3**, **OA1-T3** (count factor²),
**CD1-D1** ((in-out/round)² = (4 − d) − d·m), **RC1** (rungs; r/r_s = 3q²/(q² − p²)), **PT1** (exact turns),
**SW2** (turning centre), **QP1** (three rates), **TH1** (rung 1/3).

## Results

**F1 — every circle, d cuts, centre at rest.** Count factor² = 1 − d·m/2, and

```text
(in-out rate / round rate)²  −  (seen − lost)  =  3 − d .
```

In three cuts, and only there, **F = seen − lost = (in-out / round)²** on every free circle. AC1-A3 is the case
F = 0.

**F2 — a right triangle.** Since seen + lost = 1, in three cuts

```text
seen²  =  lost²  +  (in-out / round)² .
```

On the rung p/q (RC1) the three sides are (q² − p², 2pq, q² + p²)/2q²: a whole-number right triangle, the same
triples as PT1's exact turns. The rung 1/3 (TH1, QP1) is the triangle 3, 4, 5: seen 5/9, lost 4/9, flip 1/9.

**F3 — a turning centre.** On SW2's form, with y = a/r^(3/2) (lengths in r_s/2):

```text
seen = (1 − 3/r + 2y) / (1 + y)²
½[(in-out)² + (up-down)²] / round²  =  seen / (1 − a·Ω)²
```

At rest the up-down rate equals the round rate and this is F1. With a turn the last stable circle is no longer at
seen = lost:

```text
turn a        −1         −0.5     0      0.29     0.5      0.9      0.99
seen          108/169    0.582    1/2    0.431    0.363    0.138    0.027
```

near rest, seen = ½ − (√6/12)·a. The diagonal statement of AC1-A3 holds for a centre at rest only.

**F4 — a measured circle.** For GRO J1655−40 (441 ± 2, 298 ± 4, 17.3 ± 0.1 Hz) the three rates alone give
½[(in-out)² + (up-down)²]/round² = 0.5141 ± 0.0031, and with QP1's mass, turn and radius

```text
seen = 0.4931 ± 0.0030 ,     (in-out / round)² = 0.105 ± 0.006 .
```

The circle sits 2.3 errors from seen = lost. The line gives no reason for it to be there exactly: at rest the rung
1/3 has seen 5/9, and the turn of this centre moves it. Recorded as a number, not as a result.

## What is put in

- QP1's rates and SW2's form for F3; the measured rates of F4 as read in QP1.

## What is not shown

- F1 and F2 are short algebra on CD1 and OA1; their content is that the law's F is a measured ratio, and that the
  identity singles out three cuts (the same fact as CD1's closure and AC1-A3, now on every circle).
- F3's formulas follow from the known rates of a turning centre; not claimed new.
- A form of "seen" and "lost" for the turning centre in which the last stable circle is again on the diagonal:
  not found.

## Claim boundary

```text
(IN-OUT / ROUND)² − (SEEN − LOST) = 3 − d                                            PROVED (d = 3 … 7)
THREE CUTS: SEEN² = LOST² + (IN-OUT / ROUND)² ; RUNGS ARE WHOLE-NUMBER RIGHT TRIANGLES  PROVED
TURNING CENTRE: MEAN CROSS RATE = SEEN / (1 − aΩ)² ; EDGE AT ½ − (√6/12)a + …        PROVED
LAST STABLE CIRCLE ON THE DIAGONAL FOR A TURNING CENTRE                              REFUSED (exact table)
```

## Reproduce

```text
python fk1_flip_of_a_free_circle.py
python -m unittest test_fk1
```
