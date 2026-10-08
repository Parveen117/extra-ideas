# EN1 — the native mass bound on wave energy holds in all three sectors; for a boost it is the Doppler square

Stage of the physics line. Exact rationals, stdlib only: `en1_native_energy_bound.py`, `test_en1_exact.py` (6 tests).

## Sources used (read, unchanged)

| Source | Statement used |
|---|---|
| RH-Framework `T01_NATIVE_SUMMABILITY_AND_RECOGNITION_COMPLETION.md` (PROVED) | gauge D(r + tι) = \|r\| + \|t\| (2.1); T01-A M(a⋆b) ≤ M(a)M(b) (2.4); energy E = Σ N (4.3); weighted cut-square identity (5.1); T01-C E(a⋆ξ) ≤ M(a)²E(ξ) |
| this line | WQ1-W1 (a wave mode is a pair, energy (ω/2)(w² + p²)); EM1; GW1-V1; RB1; LT1 |

The ledger proves T01 over the cut field with ι² = −1. The line works with three sectors, u² = −1, 0, +1.
The line's wave energy (WQ1) is the ledger's recognition energy of the pair: (ω/2)·N(w + pι).

## Statements

Write a scalar as z = r + tu, D(z) = |r| + |t|, own square N_u(z) = z†z = r² − u²t², and
E��(z) = r² + t² for the positive energy of the two readings.

**T1 (T01-A in every sector).** D(zw) ≤ D(z)D(w) and M(a⋆b) ≤ M(a)M(b) for u² = −1, 0, +1.
*Proof.* zw = (rs + u²·tv) + (rv + ts)u; the ledger's proof of (2.3) uses only |u²| ≤ 1. ∎

**T2 (the own square bounds nothing outside the turn sector).** N_u is multiplicative in every sector,
but for u² = +1 it is not positive: 1 + u and 1 − u have square 0 and their sum has square 4.
No bound E(a⋆ξ) ≤ c·E(ξ) can hold with the boost sector's own square.

**T3 (T01-C in every sector, with the positive energy).** For finite path elements over any of the three sectors,

  E₊(a⋆ξ) ≤ M(a)² E₊(ξ).

*Proof.* For scalars, |rs + u²tv| ≤ |r||s| + |t||v|, and identity (5.1) with weights (|r|, |t|) gives
(|r||s| + |t||v|)² ≤ D(z)(|r|s² + |t|v²); likewise (|r||v| + |t||s|)² ≤ D(z)(|r|v² + |t|s²).
Adding, E₊(zw) ≤ D(z)²E₊(w). For paths the ledger's proof of T01-C applies unchanged, since it uses only
this scalar bound and (5.1) with the weights D(a_α). ∎

**T4 (boost: the bound is the Doppler square, reached only by a null reading).** For the unit boost
a = (1 + βu)/N, N² = 1 − β²:

  M(a)² = (1 + β)/(1 − β),

  1/M(a)² ≤ E₊(aw)/E₊(w) ≤ M(a)²,

with the upper value exactly on the co-moving null reading w ∝ 1 + u, the lower exactly on 1 − u,
and strict inequality for every other reading.

**T5 (turn and shear).** A unit turn keeps E₊ exactly (gain 1) although its mass bound is larger
((3 + 4ι)/5: bound 49/25). A shear 1 + βu, u² = 0, never reaches its bound (1 + |β|)².

**T6 (along the tower).** With β′ = 2β/(1 + β²) (LT1), M(a′)² = (M(a)²)². From r = 9 r_s the bound
for the gravity element is 2, 4, 16, 256, … : it squares at each level and is unbounded at the horizon.

## Reading

The ledger's constant M(a)² is not slack in the boost sector. It is the factor by which a boost can
multiply the energy of a wave pair, and only light going the same way gets all of it. EM1 had this
factor as a change of frame; here it is the mass of the frame element, from the RH ledger's own bound.
GW1-V1's "difference of two squares" is the boost sector's own square N_u of (rate, gradient); T2 says
why it cannot serve as an energy and T3 which square can.

## What is put in, what is not claimed

* T1 and T3 extend PROVED ledger statements to u² = 0, +1 by the ledger's own proofs. The completions
  (Cauchy presentations, T01-A/B) are not re-derived for those sectors here; only the finite inequalities.
* Random exact trials in the code are checks, not the proof.
* No new measured number: the Doppler factor and the divergence at the horizon are known.
* E₊ is a choice of positive square for a pair of readings in a boost-sector element; it depends on the
  frame, as energy does. Nothing frame-independent is claimed for it.
