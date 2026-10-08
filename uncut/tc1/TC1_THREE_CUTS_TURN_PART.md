# TC1 — with three cuts, what the self-dagger first cut leaves out is a second four-reading

Stage of the uncut line. No measurement. Exact (sympy): `tc1_three_cuts_turn_part.py`, `test_tc1.py` (5 tests).

Continues FC1, which found on the two-cut carrier one missed quantity W with det = m² + W². The question left
open there: where that part sits on the three-cut carrier, where ι commutes with the cuts.

## Sources used (read, unchanged)

extra-ideas `physics/fr1` (three cuts C₁ = K, C₂ = RK, C₃ = ιR; no fourth), `in1` (reading tensor; det = n² − r·r),
`em1` (F = E + ιB; F·F = (E·E − B·B) + 2ι E·B), `dm1` (the third cut is the unseen leg), `gb1` B1 (a reading is blind to the
phase turn), `ob1` (content counted along the phase); `uncut/fc1`, `up7`.

## Set-up

Every element of the three-cut carrier is, uniquely,

  ρ = h + ι k,  h = n + r·C,  k = d + e·C,

with h and k each their own dagger. A reading in the sense of IN1 is the case k = 0. A change of frame is
ρ → GρG†, det G = 1. Write a·b = a₀b₀ − a⃗·b⃗.

## Statements

**T1 (two four-readings).** Under every change of frame h and k transform separately, each as IN1's reading.
There are three invariants: h·h, k·k and h·k. (Two cuts had two: m² and W.)

**T2 (the determinant).** det ρ = (h·h − k·k) + 2ι h·k — the same form as F·F = (E·E − B·B) + 2ι E·B for the
light field. The light field is the traceless case, with E as h and B as k: for light, the part a self-dagger
cut leaves out is the magnetic field.

**T3 (the two-cut case inside it).** The turn of the two-cut carrier is R = −ιC₃. FC1's W is minus the third
component of k, and det = m² + W² is T2 with k along C₃. The missed part points along the cut that the two-cut
frame does not have: for that frame it is a number that no change of its frame alters (FC1-T1); for the
three-cut frame it is one component of a direction.

**T4 (the time part of k is a phase).** A change of frame by a phase does nothing to any element (GB1-B1).
Multiplying the element by Exp(ιχ) turns h and k into each other through χ. So the time part d of k is not
a direction; it is the place of the element on the phase turn, the direction along which OB1 counts content.

**T5 (at rest).** With h = (m; 0): h·k = m·d. When d = 0 the determinant is real and

  M² = m² + |e|²,

with e a direction in space that turns with the frame; for a moving element k gains a time part and h·k stays 0.
Opposite e cancel in a pair, and the pair is a plain reading.

## Answer to "did a dimension not exist then"

Yes, in this sense. On two cuts the missed part is the single number W. On three cuts it is a four-reading k,
and W was its component along the third cut — a direction the two-cut frame has no cut for. Its space part is
a direction carried by the element at rest, adding from system to system and cancelling in opposite pairs; its
time part is a phase. Neither is seen by an energy reading.

## What is put in, what is not claimed

* The reading of h as energy and momentum is IN1's. No name from physics is given here to k: that its space
  part behaves as a direction at rest orthogonal to h, and its time part as a phase, is what is shown. Whether
  these are the spin and the charge of a body is not derived, and nothing is counted in whole units.
* T2's identification for light uses EM1's F = E + ιB.
* No physical system, no number.

## Open gates

1. Counting: QC1's rule (a content is silent iff its return is a whole turn) applied to k — whether e and d come
   in whole units, and of what.
2. The loop reading of e: the return around a loop for an element with k ≠ 0 on three cuts (RV1's gyroscope).
3. The relation M² = m² + |e|² for composite elements: how |e| and m change when elements combine (MS1).
