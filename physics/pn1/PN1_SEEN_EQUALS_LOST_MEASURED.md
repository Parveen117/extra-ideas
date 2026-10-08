# PN1 — Pure numbers of the centre at rest, against measurement: seen = lost holds to two parts in 100,000

Monty Dabas. 9 October 2026. Python 3.12. Symbolic algebra (sympy) and plain numbers.

DO1: the observer on the diagonal loses what he sees. This stage collects the pure numbers that the line's field
of a centre at rest gives — ratios with no unit and no constant in them — and puts them against measured values.

Sources read before building: **MO1** (the form; M4 light circles), **MA1** (N·S = 1), **LT1** (r = s(1 + r_s/4s)²),
**LN1-N5 / GB1** (seen = lost in the bending of light), **DO1**, **MC1-T5** (spreading the clock factor instead of
the memory differs at second order), **CD1-D2**, **DG1**. Measured values: **arXiv:1409.7871** (text and Table IV)
and **arXiv:2311.08680**, both read at source.

## Results

**W1, W2 — the two weak-field numbers.** In the even layout of LT1 the time part of the form is
((1 − x′)/(1 + x′))² and the space part (1 + x′)⁴, x′ = r_s/4s. With U = r_s/2s:

```text
time part    1 − 2U + 2·b·U² + …        b = 1
space part   1 + 2·g·U + …              g = 1 :   lost / seen = 1
```

g is the ratio of what the space part and the clock part each give to the bending of light: LN1-N5's
"seen = lost", DO1's diagonal. b = 1 says the memory is exactly r_s/r.

**W3 — the other choice.** Had the clock factor been spread at least cost instead of the memory (MC1-T5), the
orbit ratio would be 1 − (5/2)·x: five sixths of the advance, 35.8″ per century for Mercury.

**W4 — light.** The light circle is at (3/2)·r_s, its reach is (3√3/2)·r_s, and the dark disc seen from far has
diameter 2√27 × GM/c²D = 10.39 × GM/c²D.

**W5.** At the last stable circle the binding is 1 − √(8/9) = 0.0572 of the rest count.

## Measured

```text
pure number                              line        measured                          how
lost / seen − 1                          0           (2.1 ± 2.3) × 10⁻⁵                radio delay past the Sun
(1 + lost/seen) / 2                      1           0.99992 ± 0.00023                 radio sources past the Sun
second-order number − 1                  0           (−4.1 ± 7.8) × 10⁻⁵               Mercury's perihelion
4·b − g − 3                              0           (4.4 ± 4.5) × 10⁻⁴                Moon and Earth falling to the Sun
dark disc / (GM/c²D)                     10.39       9.50 ± 1.37                       centre of the Milky Way
advance of Mercury, ″ per century        42.98       42.98   (other choice: 35.8)      MO1
```

Seen = lost is measured to two parts in 100,000. The second-order number separates "memory is what is spread"
from "the clock factor is what is spread" by 40% against a measured uncertainty of 0.008%: MC1's assumption on
this point is decided by measurement.

## What is put in

- MO1's form (MA1) and the dictionary r_s ↔ GM/c². For the dark disc: the mass-to-distance ratio from the stars'
  orbits, as quoted in the source.

## What is not shown

- These numbers are the known ones for this field. The line's part: each is stated as a ratio of its own
  quantities (lost : seen; memory spread at least cost), and W3 turns an assumption of MC1 into a measured fact.
- The centre of the Milky Way turns; the disc number for a turning centre (SW2) is not computed here. The source
  gives −8% to 0 for that range.

## Claim boundary

```text
LOST/SEEN = 1 AND SECOND-ORDER NUMBER = 1 FROM THE FORM                         PROVED (symbolic)
OTHER CHOICE OF MC1-T5 ⇒ 5/6 OF THE ADVANCE                                     PROVED ; REFUSED BY MEASUREMENT
LIGHT CIRCLE REACH (3√3/2) r_s ; DISC 2√27 GM/c²D                               PROVED
FIVE MEASURED PURE NUMBERS                                                      ALL WITHIN THEIR UNCERTAINTY
A NUMBER THAT DIFFERS FROM THE KNOWN THEORY                                     NONE
```

## Reproduce

```text
python pn1_seen_equals_lost_measured.py
python -m unittest test_pn1
```
