# OA1 — The once-around turn with the full form of MO1: 2π(1 − count factor), first order (3/2)·π·r_s/r

Monty Dabas. 8 October 2026. Python 3.12. Symbolic algebra (sympy).

Third item of "Next" in the gravity thesis: the once-around turn with the full form of MO1 instead of the static
frames of GR2-G4. The thesis recorded that GR2-G4's turn, π·r_s/r, is two thirds of the measured effect on an
orbiting gyroscope and that "the remaining part is not in this model".

Sources read before building: **GR2-G4** (a frame boosted radially by the fall rapidity, carried once around,
returns turned by 2π(1/N − 1)), **MO1** (the form; M4: on a free circle (dφ/dt)² = r_s/2r³ and count factor² =
1 − (3/2)·r_s/r; light circles at r = (3/2)·r_s), **MA1** (the form follows from the law for one centre at rest),
**CV1** (the connection that keeps the form and has no torsion; T1: frame route and coordinate route agree),
**GR1** (N² = 1 − m, static acceleration m/(2rN) in units c = 1), **CL2** (free readings take the largest count).

## Setting

The count form of MO1 in the plane of the circle, dτ² = dt² − (dr + β dt)² − r² dφ², β² = r_s/r. A direction is
carried along a circle of radius r gone round at the rate Ω = dφ/dt, with CV1's connection. Two pure numbers:

```text
x = r_s / r   (the memory m at the circle) ,        y = (r·Ω)² .
```

On a free circle the carried direction is moved by the connection alone. On a held circle one rule is added: the
direction is kept across the history and is not turned otherwise.

## Results

**T0.** Along the circle the carried direction makes a pure turn (the transport has one zero rate and one pair ±iω).

**T1 — every circle.** Against the basis that goes round with the circle, the direction turns per circuit by
2π·(1 − 3x/2)/√(1 − x − y). So against far directions it turns, per circuit, by

```text
turn / 2π  =  1 − (1 − 3x/2) / √(1 − x − y) .
```

**T2 — in the line's quantities.** With N² = 1 − x, static acceleration g = x/(2rN) and local speed v (y = v²N²):

```text
turn / 2π  =  1 − γ·( N − r·g ) ,        γ = 1/√(1 − v²) .
```

**T3 — the free circle.** With MO1's period law (y = x/2):

```text
turn / 2π  =  1 − √(1 − 3x/2)  =  1 − (count factor of MO1-M4) .
```

Per circuit the carried direction falls short of the far directions by exactly the share of count that the
circling reading falls short of the far clock. To second order: (3/2)·π·x + (9/16)·π·x².

**T4 — against GR2-G4.** GR2-G4 gives 2π(1/N − 1) = π·x + (3/4)·π·x². The first order of the full form is 3/2
of it, exactly. The missing half is there once the direction is carried along the moving history of MO1 instead
of being passed between static frames.

**T5 — carried slowly on a held circle (y → 0).** turn/2π = (1 − N) + m/(2N): the cone part of the thesis plus
radius × static acceleration.

**T6 — no field (x = 0).** turn/2π = 1 − γ: the turn of a direction carried round a circle at speed v, with
nothing put in for it.

**T7 — the light circle.** At x = 2/3 (r = (3/2)·r_s, MO1-M4's light circles) the turn is one whole turn per
circuit for every speed: the carried direction keeps step with the circle.

## Numbers

```text
circuit of the Earth at r = 7020 km (the thesis's example)
    GR2-G4                       0.8188 milliarcsecond per circuit
    full form (T3)               1.2282 milliarcsecond per circuit   →  6.621 arcsec per year
    measured (Gravity Probe B)   6.6018 ± 0.0183 arcsec per year      arXiv:1105.3456, abstract, read at source
circuit of the Sun at 1 astronomical unit (the Earth–Moon pair as the carried direction)
    full form (T3)               19.19 milliarcsecond per year
```

The formula at 7020 km is 0.29 % above the measured value (1.05 of its stated uncertainty). The rate goes as
r^(−5/2), so the measured band corresponds to r = 7028 ± 8 km. The paper gives the altitude of the orbit, 642 km,
and not its mean radius; a semi-major axis of 7027.4 km is recalled for the mission (from memory, not verified
at source) and would give 6.604. The Earth's flattening and the Sun's field are not included here. For the
circuit of the Sun, the measured value recalled from lunar ranging is about 19 milliarcsecond per year (from
memory, not checked at source).

## What is put in

- The form of MO1 (derived in MA1 for one centre at rest) and CV1's connection.
- For the free circle: nothing else. The period law is MO1-M4.
- For held circles (T2, T5, T6): the rule that the direction is kept across the history and not turned otherwise.
- For the numbers: GM of the Earth and of the Sun, the radius of the circuit, the length of the year.

## What is not shown

- This is the known turn of a gyroscope on a circular orbit (geodetic turn), reached with the line's form and
  connection. The second-order coefficient 9/16 is the known one for a circle; it is far below measurement.
- Circles only. Other orbits, and a centre that turns (see SW1), are not in this stage.
- T3's reading "shortfall of turn = shortfall of count" is an exact identity on free circles of this form; no
  general law is claimed from it.

## Claim boundary

```text
TURN PER CIRCUIT ON ANY CIRCLE: 1 − (1 − 3x/2)/√(1 − x − y)               PROVED (symbolic)
FREE CIRCLE: 2π(1 − √(1 − 3x/2)) = 2π(1 − COUNT FACTOR OF MO1-M4)          PROVED (symbolic)
FIRST ORDER = 3/2 × GR2-G4 ; SECOND ORDER 9/16·π·x²                        PROVED (symbolic)
THESIS: "THE REMAINING PART IS NOT IN THIS MODEL"                          ANSWERED: it is in MO1's form
AGREEMENT WITH THE MEASURED GYROSCOPE TURN                                 0.3 % at r = 7020 km; orbit radius and flattening not treated
ANYTHING NEW PREDICTED                                                     NO
```

## Reproduce

```text
python oa1_once_around_with_the_full_form.py
python -m unittest test_oa1
```
