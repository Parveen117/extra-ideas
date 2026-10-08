# TT1 — where time sits in the thermo diagram

Stage of the uncut line. No measurement and no physical constant. Exact (sympy):
`tt1_time_in_the_thermo_diagram.py`, `test_tt1.py` (6 tests).

The owner's statement: relativity has become a special case, connected even before thermo. Now use the thermo
space together with time in place of space, and see what relation between time and thermo comes out.

## Sources used (read, unchanged)

| Source | Statement used |
|---|---|
| response-geometry RMG10 T2, T4, T6 | the response space is SL(2,R) with a metric of signature (2,1); u, v spacelike, t = w/m timelike; causal sectors; the tower keeps t/ρ |
| response-geometry RMG11 T3, RMG9 T1–T2 | the three native exponentials; L = m(1 + uK + vS) + wR |
| extra-ideas `uncut/up4`, `up6`, `up7` | w; the cross-corner pair; cuts u² + v² − w² = 1 and the four-step walk |
| extra-ideas `physics/rb1`, `lt1`, `gb1`, `tl1` | return = boost; tower generation x → 2x/(1 + x²); unit block Exp(ψn); units follow the clock |

## Statements

**T1 (the space of cuts has two space directions and one time direction).** For a cut κ = uK + vS + wι,
u² + v² − w² is kept by every change of frame on the carrier. Exp(θι) turns (u, v) by 2θ and leaves w.
Exp(ηK) and Exp(ηS) mix w with v and with u as boosts of rapidity 2η. The turn-part w of the cut is the
time direction; the self-dagger cuts — the closed diagrams of UP7 — are the instant w = 0.

**T2 (one cycle of the diagram is a boost).** Write w = sinh η. Then

  (ικ)² = Exp(−2η n),

with n a self-dagger cut. Going once round the diagram — κ, ι, κ, ι — is the identity when w = 0 and
otherwise a boost of rapidity 2η. Cycles add: k cycles are Exp(−2kη n). In the reading of RB1 and GB1
(a unit block Exp(ψn) has speed tanh ψ and clock factor 1/cosh ψ):

  speed of one cycle = tanh 2η = tower(tanh η),  clock factor = 1/cosh 2η = 1/(1 + 2w²).

For w = 3/4: tanh η = 3/5, one cycle has speed 15/17 and clock factor 8/17. One cycle of the diagram is one
generation of the λ-tower (LT1's squaring) applied to tanh η.

**T3 (the thermo value).** For a response element L = [[A, B₁], [B₂, C]],

  tanh η = (B₂ − B₁) / √((A − C)² + (B₁ + B₂)²) = t/ρ :

the antisymmetric part of the response over the anisotropy of its symmetric part. Three cases, RMG10's
causal sectors:

| ((A − C)/2)² + B₁B₂ | the response's cut | its cycle |
|---|---|---|
| > 0 | a cut with turn-part w | a boost of rapidity 2η |
| = 0 | degenerate | — |
| < 0 | none: the traceless part is itself a turn | a rotation |

A reciprocal response (B₁ = B₂) has η = 0: no time part.

**T4 (the tower does not change it).** tanh η is the same at every level L → L + λL².

**T5 (an example with nothing put in by hand).** For a gas with constant capacities and the cross-corner pair
of UP6 (one cut at fixed V, the other at fixed T): A/U = S(S + c)/c², B₁/U = R/c, B₂ = C = 0, and

  tanh² η = (Rc)² / (S²(S + c)² + (Rc)²).

Large S: η → 0. S → 0: tanh η → 1, the light cone.

## Reading

In the thermo diagram, time is not an extra axis laid beside S and V. It is already there as the third
coordinate of the cut: the part of the cut that is a turn. Where the diagram closes (Maxwell's one equation,
w = 0) there is no time part. Where it does not close, one trip round it is a boost, and trips add like
rapidities — a count of cycles is a clock. The relation asked for is T3: the rapidity per cycle is fixed by the
ratio of the non-reciprocal part of the response to its anisotropy.

## What is put in, what is not claimed

* "Speed" and "clock factor" in T2 are the readings of a unit block used in RB1 and GB1. That a cycle of a
  thermodynamic process is a boost of a physical clock is that dictionary, not a derivation. No second, no
  metre and no constant appears.
* T5 depends on the zero from which S is counted: the scale operation S·∂/∂S is not unchanged by a shift of S.
  The same is true of λ_p = −ST/C_p. The example shows the form of the relation, not a number for a gas.
* TL1's statement (units follow the clock) would turn the clock factor 1/(1 + 2w²) into a factor on the local
  unit of fluctuation. That step is not taken here.
* A process in time — a state moving on the thermo plane — is the river of RV1 on that plane, with the
  non-reciprocal part as its spin (RV1-T9). It is not rebuilt here.

## Open gates

1. The process form: the frame e₀ = ∂_t + v·∂ on the thermo plane with v = −L X, its count, and whether its
   clock factor is T2's.
2. The zero of S: the relation of T3 for a pair of cuts that does not depend on it.
3. The sector D < 0: a response with no cut — what replaces the diagram there.
