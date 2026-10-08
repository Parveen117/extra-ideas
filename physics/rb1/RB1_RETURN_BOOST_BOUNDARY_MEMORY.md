# RB1 — the recognised return is a boost; the horizon end keeps boundary memory, the flat end does not

Stage of the physics line. Exact rationals, stdlib only: `rb1_return_boost_boundary_memory.py`, `test_rb1_exact.py` (8 tests).

## Sources used (read, unchanged)

| Source | Statement used |
|---|---|
| Recognition-Kernel-Framework `research/recognition_return` R1 | RR1 first-return law F = 1 + C F B F; RA1–RA3 completed return F = 1 + dL, d = q(1 − d²); RC1 strong return selects (1 + L)/2 |
| same, `r2` R2 | RD1 pairing = closure predicate; RD2 backward recursion x_j = a_j/(1 + a_j x_{j+1}); RI1 exact interval; RI2 tail forgotten iff Σ 1/a_j = ∞ |
| Publications `papers/native-critical-response` (research branch) | CR-1 channel split of the inverse, CR-2 source/observer selection rule |
| this line | GB1 gravity element (1 + β r̂·C)/N, LT1 tower generation x → 2x/(1 + x²), MO1 fall speed β² = r_s/r |

None of these return theorems had been used in the physics line before this stage.

## Statements

**T1 (return = boost).** For the uniform chain of R1 with coupling q, the completed return is
F = 1 + xL with q = x/(1 − x²). Writing x = tanh η: F = Exp(ηL)/cosh η and **2q = sinh 2η**.
So twice the coupling is the proper speed x′/√(1 − x′²) of the *next* tower level x′ = 2x/(1 + x²) (LT1).
Returns compose by the speed-addition law: (1 + xL)(1 + yL) = (1 + xy)(1 + (x + y)/(1 + xy) L).
Benchmark: q = 6/5 ⇒ x = 2/3, unit norm 5/9, next level 12/13.

**T2 (the gravity element is a chain return).** Along the tower x_{j+1} = 2x_j/(1 + x_j²) the cells
a_j = x_j/(1 − x_j x_{j+1}) are positive, the chain closed by its own next value returns x_0 exactly, and

  1/a_j = (1 − x_j²)/(x_j (1 + x_j²)),  1 − x_{j+1}² = ((1 − x_j²)/(1 + x_j²))².

The clock factor N² = 1 − x² squares (or better) at each level.

**T3 (the horizon end remembers).** Hence Σ 1/a_j is finite going inward, and by R2 RI2 the interval of
possible outer returns does not close. From x_0 = 1/3 (r = 9 r_s): Σ 1/a_j < 4, and the outer return lies
anywhere in [0.27305, 0.39571] (width 0.122665, stable after five levels) depending on what closes the
horizon end. Closure by the null cut 2P₊ (tail 1, the cut RC1 selects) gives 1/3 — the point-mass profile;
a wall (tail 0) gives one end of the interval,
0.39571 for an odd number of cells and 0.27305 for an even number (R2 RI2's two limits).

**T4 (the flat end forgets).** Going outward with x_j → 0, every 1/a_j ≥ 1/x_j − 1, so Σ 1/a_j = ∞ and
the reading is independent of what closes the far end (shells r = (j + 3)² r_s: interval width below
10⁻¹³⁰ after 50 cells).

**T5 (the pole at closure is seen by one channel only).** (1 + xL)⁻¹ = P₊/(1 + x) + P₋/(1 − x).
At x = 1 the return is 2P₊; the inverse has a simple pole carried by P₋ alone; the P₊ reading stays 1/2;
R exchanges the channels (R P₊ = P₋ R). With x² = r_s/r the pole factor 1/(1 − x²) is the radial factor
of the static reading: singular for the observer who retains P₋, regular for the pure P₊ (infalling) reading.

## What is put in, what is not claimed

* The identification "x = fall speed β" is a dictionary between two unit elements of the same algebra
  (R1's 1 + xL and GB1's 1 + β r̂·C). It is not derived from a source law.
* The shell sequence in T3 is the tower of LT1. Any inward sequence that reaches x → 1 in infinitely many
  cells with summable 1 − x_j² gives the same verdict; a different spacing with that summability changes the numbers, not the verdict.
  Later note (HM1-T4): a spacing without it — for example x_j = 1 − 1/(j + 2) — reverses the verdict.
* No measured number is predicted. For a static field a different horizon closure is read outside as a
  different r_s. Whether a time-dependent probe can separate closures is not addressed here.
* The inverse of a return is not a propagator (CR-1's own caveat is kept).

## Open gate named by this stage

Which closure does a physical horizon have? RC1 selects 2P₊ only in the joint limit (coupling → ∞ with
depth/coupling → ∞). A closure other than 2P₊ would be the first place where the framework's exterior
reading can differ from the point-mass one with the same cells.
