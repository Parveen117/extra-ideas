# BC1 — The boundary and the count: where whole numbers come from

Monty Dabas. 8 October 2026. Python 3.12, standard library only. Exact
rational arithmetic.

LT1 found that the λ-tower ends on the boundary where powers add
exactly, and left open whether that is the counting of waves (WQ1).
This stage checks.

Sources read before building: **RMG2-T1(c)** ((RH)² = −Δ: complex
inside the disc, dual on the boundary) and its holonomy
Hol = Exp(Θ J_H), J_H = RH/√Δ; **LT1**; **WQ1**; **QC1-T2/T3**;
**TL1-L1**; **GR1**; **SY1**.

## Results

**B1 — inside: a return is a turn.** (RH)² = −Δ. With the
un-normalised generator,

```text
Exp( s·RH ) = cos(s d)·1 + ( sin(s d)/d )·RH ,       d = √Δ .
```

A turn through the angle s·d. Returns compose as turns; a whole turn
needs s = 2π/d. Whole numbers (QC1's silence, WQ1's whole area) need
this period.

**B2 — on the boundary: a return is an addition.** At Δ = 0,

```text
Exp( s·RH ) = 1 + s·RH ,       (1 + s₁N)(1 + s₂N) = 1 + (s₁ + s₂)N ,
```

and no s other than 0 comes back to 1. There is no whole turn there:
no silent content, no step. The counting is exact but it is plain
addition.

**B3 — the approach.** For the gravity element at clock factor N the
trace-normalised invariant is Δ_n = N², so d = N: the turn per unit of
far parameter is the clock factor (TL1-L1). Seen from far, the unit of
a whole turn is 2π/N and the step of a wave's energy is κωN. The count
n stays a whole number in every frame (WQ1-W4); the size of a step,
seen from far, goes to zero at the horizon.

**B4 — the bracket is additive too.** For a pair, the loop

```text
(1 + aE₁₂)(1 + bE₂₃)(1 − aE₁₂)(1 − bE₂₃) = 1 + ab·E₁₃
```

returns the enclosed area in a central slot that only adds. The area of
a wave is a quantity of the same additive kind as the boundary's.

## Answer

They are not the same thing; they are the two halves of it.

```text
whole numbers  =  an additive quantity (dual sector: B2, B4)  read as a turn (circular sector: B1)
```

The boundary supplies the addition, the interior supplies the period.
On the boundary alone the period is infinite and nothing is stepped.
So the end of the tower is not where quantisation hides; it is where,
seen from far, the steps close up and counting becomes continuous
addition. Quantisation needs both sectors at once — a pair (additive
area) inside the disc (a turn).

## What is not shown

- B3 uses the clock factor as the scale of the turn; this is the known
  slowing of clocks, in the tower's terms. Nothing new is predicted.
- "Steps close up at the horizon" is a statement about energies carried
  to far away. Locally the step is κω as everywhere.
- A count belonging to the horizon itself is not derived, and B2 gives
  no support for one: the boundary sector has no period of its own.

## Claim boundary

```text
(RH)² = −Δ ; RETURN = TURN THROUGH s√Δ INSIDE THE DISC                    PROVED
RETURN = ADDITION ON THE BOUNDARY; NEVER COMES BACK                       PROVED
Δ_n = N² FOR THE GRAVITY ELEMENT; FAR STEP κωN                            PROVED
LOOP OF A PAIR RETURNS ITS AREA ADDITIVELY                                PROVED
WHOLE NUMBERS = ADDITIVE QUANTITY READ AS A TURN                          FOLLOWS from the four
QUANTISATION LOCATED AT THE BOUNDARY ALONE                                REFUSED (B2)
A COUNT OF THE HORIZON                                                    NOT DERIVED
```

## Reproduce

```text
python bc1_boundary_and_count.py
python -m unittest test_bc1_exact
```
