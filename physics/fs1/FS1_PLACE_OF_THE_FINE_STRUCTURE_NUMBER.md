# FS1 — Where 1/137 sits in the line, exactly; and what does not give it

Monty Dabas. 9 October 2026. Python 3.12. Symbolic algebra (sympy) and plain numbers.

Every pure number the line has given so far is a ratio of counts (ONE_LAW.md). The synthesis lists one pure
number of another kind as open: the size of one unit of light's content, about 1/137. This stage does not derive
it. It finds its exact place and records what fails.

Sources read before building: **GM1** (charged displaced centre: memory Re(r_s/R) − q₂/R·R̄; moment/charge = a;
ratio 2; shortfall 0.00116), **SW2** (swirl constant r_s·a), **EG1** (q₂), **EM1** (the wave gives the place of the
unit; the number to aim at is the ≈ 1/137 coupling), **OB1** (content q), **RC1** (closure and the turn left over),
**DO1**, **SC1** (no law in the block gives a constant with a unit).

## Results

**A1 — the one number of the centre without the unit of gravity.** The centre has three constants: r_s (mass),
a (displacement), q₂ (charge). In

```text
q₂ / ( r_s · a )  =  k_e e² / ( 2·J·c )
```

G and the mass cancel. For turning content J = ħ/2 it is k_e e²/ħc: the fine-structure number. With the electron's
constants: 1/137.036.

**A2 — in the far field.** Memory plus ι × swirl potential of the charged displaced centre, on its axis:

```text
r_s / d   +   r_s·a·( ι − α ) / d²   +  …
```

At second order the turning enters along ι and the charge along −1, in the ratio 1 : α. The number 1/137 is the
size of the charge term measured against the turn term of one and the same centre.

**A3 — the left-over turn.** GM1's ratio 2 means the moment's direction keeps step with the circuit in a magnetic
field: closure, in RC1's sense. The measured ratio leaves 2π·(g/2 − 1) = 0.007286 radian per circuit, which is α
to within 0.15%. To first order the fine-structure number is the turn left over per circuit, in radians.
(The relation between the moment ratio and the turn per circuit is recalled, not derived in the line.)

**A4 — what does not give it.**

```text
the diagonal on this centre (charge term = turn term)        would be α = 1                        refused by measurement
small ratios n/m ≤ 12 with powers of π and √2                closest is 0.4% off                   refused
the known near miss 4π³ + π² + π                             2.2 × 10⁻⁶ off; measured to 10⁻¹⁰       refused
```

## What this settles

```text
EM1: "the place of the unit"        the place is exact: α = q₂ / (r_s·a) for a centre of turning ħ/2
the kind of number                  not a ratio of counts of cuts or planes, in any of the simple forms tried
what a derivation must be           a law that fixes the charge term of a centre against its turn term
```

## What is not shown

- A1–A2 are identities of the dictionary; they do not fix the value.
- The left-over turn of A3 is not computed from the law. GM1's centre gives exactly 2.
- A4 is a finite search; it does not prove that no expression exists.

## Claim boundary

```text
q₂/(r_s a) = k_e e²/(2Jc) ; = α FOR J = ħ/2                                      PROVED (symbolic) ; 1/137.036 from the constants
FAR FIELD: SECOND COEFFICIENT r_s a (ι − α)                                      PROVED (symbolic)
MEASURED LEFT-OVER TURN PER CIRCUIT = α TO 0.15%                                 NUMBERS (measured value recalled)
α FROM THE DIAGONAL OR FROM SMALL RATIOS OF COUNTS                               REFUSED
THE VALUE 1/137                                                                  NOT DERIVED
```

## Reproduce

```text
python fs1_place_of_the_fine_structure_number.py
python -m unittest test_fs1
```

## Later note (BR1, 9 October)

In a bound record α is the speed of the first recordable circle, the square root of its memory, and the ratio of
the charge's three lengths (BR1). Its value remains open. Nothing above is changed.
