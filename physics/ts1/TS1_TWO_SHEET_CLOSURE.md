# TS1 — a second sheet with its own flat end removes the closure freedom of the horizon

Stage of the physics line, continuing RB1. Exact rationals, stdlib only: `ts1_two_sheet_closure.py`,
`test_ts1_exact.py` (5 tests).

## Sources used (read, unchanged)

| Source | Statement used |
|---|---|
| Recognition-Kernel-Framework `research/recognition_return` R1 | RC1 strong return selects (1 + L)/2 with defect −d/(4q); RC2 cutoff and coupling limits do not commute |
| same, `r2` | RI1 interval, RI2 criterion Σ 1/a_j, SY2 cells of a profile |
| this line | RB1 (T2–T4), LT1 tower generation and radius map, PR2 two sheets |

## Statements

**T1 (the tower has a mirror).** tower(x) = tower(1/x) exactly; in radius the LT1 map (r + r_s)²/4r is
unchanged under r → r_s²/r. One tower level cannot tell the outside from its inversion: the horizon is
the fixed point of a mirror, and behind it the tower places a second copy of the outside.

**T2 (closure freedom: one sheet 0.1227, two sheets 0).** From x₀ = 1/3 (r = 9 r_s), six tower levels:

| chain | interval of possible outer readings |
|---|---|
| one sheet, any closure at the horizon end | width 0.122665 |
| mirrored second sheet, nothing behind it | 0.023195 |
| second sheet + 5 flat cells | 9.7·10⁻⁹ |
| second sheet + 20 flat cells | 5.9·10⁻⁴³ |
| second sheet + 80 flat cells | 1.2·10⁻²⁴⁵ |

The point-mass value 1/3 lies in every interval. With the second sheet's flat end, Σ 1/a_j diverges
(RB1 T4), so by RI2 the outer reading is unique and equals the point-mass one whatever closes the far end.

**T3 (the throat is R1's strong return).** The cell joining the two sheets is exactly R1's uniform coupling
q = x/(1 − x²) at the throat value, and (1 + xL)/2 misses being a cut by −x/(4q) (below 6·10⁻²⁰ at six levels).
The cut 2P₊ that RB1 had to put in as a closure is what the throat cell returns.

**T4 (RC2 seen on the chain).** The strong cell alone, with a wall behind it, returns q (unbounded), not the
cut. With the second sheet behind it, it returns the throat value. The depth that RC2 demands
(depth/coupling → ∞) is supplied by the second sheet.

## Reading

RB1 left one free datum at the horizon. It is free only on one sheet. If the element lives on two sheets
(PR2: mass and antimass sheets; FR1: with three cuts mass exists only as coupling to the conjugate sheet),
the datum is fixed and the fixed value is the point-mass profile. Conversely, a one-sheet object with a
wall (closure 0) would read 0.2730 or 0.3957 at the same cells instead of 1/3 — which of the two depends on
the parity of the number of cells (R2 RI2's two limits).

## What is put in, what is not claimed

* The second sheet's profile is taken as the mirror image of the first. T1 makes that the tower's own
  choice, but a different second-sheet profile is not excluded; any profile with a flat end gives T2's
  verdict (unique reading), only the value then follows that profile.
* x = fall speed remains RB1's dictionary. No measured number is predicted.
* Nothing here says the second sheet is traversable or physical beyond being the closure.

## Open gate

A static reading cannot separate "two sheets" from "one sheet closed by 2P₊" — the values agree to the
defect in T3. A separation needs a probe that changes the cells by a common multiplier (the calibrate-two,
predict-third protocol of the native-critical-response packet); which physical probe acts that way on the
gravity cells is not identified.
