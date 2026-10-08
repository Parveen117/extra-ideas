# SW1 — A centre that turns: the swirl a/r³ from the line's law, the turn of a carried direction, and a refusal at second order

Monty Dabas. 8 October 2026. Python 3.12. Symbolic algebra (sympy), with the functions of CV1 and TP1.

MA1 derived MO1's assumption (flat slices, one time) for one centre at rest and left the general case untreated.
This stage takes the next case: the centre turns about an axis.

Sources read before building: **MO1**, **MA1** (flat slices and one time follow from the law for a centre at
rest; frames at rest far away), **CV1** (frame curvature as order defect; empty space = zero contracted
curvature; the connection from the order defects alone), **TP1** (the law selected by equivalence of local
frames is that curvature law), **GR2-G1** (the locally pure frame falls inward with β² = m), **OA1**.
Checked that no stage of the line treats a turning centre (dock, synthesis, unused-results map).

## The family

The fall frame with one more entry: the fall history also goes round the axis, at a rate W(r) on each shell.

```text
e₀ = ∂_t − β ∂_r + W(r) ∂_φ ,   e₁ = ∂_r ,   e₂ = (1/r) ∂_θ ,   e₃ = (1/(r sin θ)) ∂_φ ,     β² = r_s/r .
```

Slices flat and one time are put in here, and then tested against the law (T3).

## Results

**T1 — the law at first order in the swirl.** The contracted curvature has, at first order, only the components
(time, round) and (radial, round), and both are the same operator:

```text
( r⁴ W′ )′ = 0      ⇒      W = W∞ + a / r³ .
```

**T2 — a constant swirl is no field.** With W constant every contracted component vanishes to all orders, and
without the centre the whole curvature vanishes: it is flat space read in a turning frame. So, as in MA1, the
frames at rest far away fix W∞ = 0, and **W = a/r³** with one constant a.

**T3 — second order: refused with flat slices.** With W = a/r³ the contracted curvature is exactly

```text
9 a² sin²θ / (2 r⁶)  ×  ( −1 ; +1 , 0 , −1 )     on (time ; radial , θ , round) ,
```

and nothing at higher order. It is not zero. So flat slices with one time hold for a turning centre to first
order in the turning and not beyond. MO1's assumption is exact for the centre at rest (MA1) and is a first-order
statement here.

**T4 — any fall velocity on flat slices.** For the frame e₀ = ∂_t + u·∇ with any u(x), and the slice directions
as the other three:

```text
the frame is in free fall                             w^i_00 = 0
along the slice the directions do not turn            w^i_jk = 0
a direction carried with the frame turns at           ½ curl u      against the slice directions, exactly
the remaining part of the connection is               the symmetric part of ∂u (stretch)
```

For the centre at rest u is radial and has no curl: a direction that falls with the frame keeps the far
directions exactly. The turn of OA1 comes entirely from the stretch met by a history that moves across the frame.

**T5 — the turning centre.** u = −β r̂ + (a/r³) ẑ × r has

```text
½ curl u  =  ( a / 2r³ ) · [ 3 (ẑ·r̂) r̂ − ẑ ] .
```

Over the poles a carried direction turns with the centre, at a/r³; in the plane of the equator against it, at
a/2r³.

**T6 — mean over a circuit.** Over a circle through the poles the mean is a/(4r³) about the axis, in the sense
of the centre. Round the axis it is −a/(2r³).

## Number

The constant a is fixed by the turning content of the centre. Here it is **put in** as a = 2GJ/c², the known
identification; the line has not derived it.

```text
the Earth, circuit over the poles at r = 7020 km
    mean turn of a carried direction        40.9 milliarcsecond per year
    seen on a direction 16.84° from the equator plane (× cos)     39.2
    measured (Gravity Probe B)              37.2 ± 7.2             arXiv:1105.3456, abstract, read at source
    value quoted there as predicted         39.2
```

J of the Earth is taken as 0.3307·M·R²·ω and the direction of the mission's guide star as 16.84° from the
equator plane; both from memory, not verified at source.

## What is put in

- The family: a swirl that depends on r only; flat slices and one time (tested in T3).
- The law of CV1/TP1 and the condition of MA1 (frames at rest far away).
- For the number: a = 2GJ/c², J of the Earth, the radius of the circuit, the direction of the gyroscope.
- For a history that moves across the frame (the orbit), the swirl part of the turn is taken from T4/T5 at the
  place of the history; the cross terms between the orbit's speed and the swirl are of relative order the square of that
  speed (10⁻⁹ here) and are not computed.

## What is not shown

- This is the known first-order field of a slowly turning centre and the known turn of a gyroscope in it,
  reached with the line's law. Nothing new is predicted.
- What replaces flat slices at second order is not found. T3 gives the exact remainder that a correct second-order
  frame has to cancel.
- A swirl that also depends on θ is not treated.
- The constant a is not derived from a source: the line has no built sources (synthesis, open list).

## Claim boundary

```text
FIRST-ORDER LAW (r⁴W′)′ = 0 ; W = a/r³ WITH FRAMES AT REST FAR AWAY        PROVED (symbolic)
CONSTANT SWIRL = NO FIELD                                                  PROVED (symbolic)
FLAT SLICES AND ONE TIME FOR A TURNING CENTRE, FIRST ORDER                 CONSISTENT WITH THE LAW
FLAT SLICES AND ONE TIME FOR A TURNING CENTRE, SECOND ORDER                REFUSED: remainder 9a² sin²θ/2r⁶ × (−1; +1, 0, −1)
CARRIED DIRECTION TURNS AT ½ curl u ON FLAT SLICES                         PROVED (symbolic, any u)
MEAN OVER A POLAR CIRCUIT a/(4r³)                                          PROVED (symbolic)
a = 2GJ/c²                                                                 PUT IN
AGREEMENT WITH THE MEASURED TURN                                           inside the stated uncertainty (20 %)
ANYTHING NEW PREDICTED                                                     NO
```

## Reproduce

```text
python sw1_turning_centre_first_order.py
python -m unittest test_sw1
```
