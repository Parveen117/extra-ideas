# QC3 — One arrow, two units: the even part is heat, the odd part is phase

Monty Dabas. 7 October 2026. Python 3.12, exact rational arithmetic only.
Continues QC2.

QC2 ended on one gate: the thermal unit sits on the real direction (κ·1),
the quantum unit on the circular one (κ·R); why? This stage answers with
the framework's own arrows. Source, read before building: **GE2-T1**
(rational arrows U_v recover the phase derivation D_v), **GE2-T2**
(typed record lift, memory M = I − S†S), **GE2-T3** (symmetric records,
ledger h = ½ tr Q, generator −L_C), all in Publications
`papers/generalized-euler-evolution`.

## 1. The split

Take one native arrow U_v and its reverse U_{−v} = U_v†. Put

```text
E = (U_v + U_{−v}) / 2        even part — the record mean of the symmetric pair (GE2-T2)
O = (U_v − U_{−v}) / 2        odd part  — first order D_v, the phase derivation (GE2-T1)
```

**T1 (hyperbolic Pythagoras, exact at every turn size).**

```text
E� − O² = I ,            M = I − E² = −O² .
```

The memory a symmetric record loses (GE2's M) is exactly minus the square
of the odd part. What the record forgets is the phase part, squared.

*Proof.* U_vU_{−v} = I and the two commute, so (E + O)(E − O) = I. ∎

**T2 (two units).** To leading order O = D_v, first order in the turn and
circular (D_v² ≤ 0); E = I + ½D_v², second order and real. Deterministic
refinement keeps O and gives the phase Exp(T·D_v); symmetric independent
records keep E and give heat with ledger τ = Σ ½|v|² (GE2-T1, T3, T4). The
circular unit is the turn; the real unit is half its square.

## 2. Circular contents: the character law

On a reading of content q (QC1) a turn θ acts as cos qθ + sin qθ·R.

**T3 (exact at finite turn).**

```text
1 − E_q = (1 − E_1) · χ_q(θ/2)² ,      O_q = O_1 · χ_q(θ) ,      χ_q(x) = sin(qx)/sin(x) ,
(1 − E_q)/(1 − E_1) ↑ q²   as the turn shrinks.
```

The heat defect of content q is the unit defect times the **square** of
the phase sum; the odd part is the unit odd part times its first power.
Under records, content q decays at rate q² per unit ledger; only q = 0
survives observation.

## 3. su(2) contents: the Casimir law, and the Yang–Mills number

With GE2's rational quaternion turns q(v) = (1 − r²/16 − v/2)/(1 + r²/16)
and the isotropic symmetric record over three axes:

**T4 (exact at finite turn).**

```text
content ½ :  record mean = χ_½/2 ,   rate (1 − S)/h = (1/4) / (1 + r²/16)      → 1/4 = (3/4)/3
content 1 :  record mean = χ_1/3 ,   rate (1 − S)/h = (2/3) / (1 + r²/16)²     → 2/3 = 2/3
ratio     :  (8/3) / (1 + r²/16)                                               → 8/3 = C_1 / C_½
```

The record mean is the normalized character; the heat rate is the
Casimir j(j+1) divided by the three axes. The numbers 3/4 and
C_1 − C_½ = 5/4 are the free reduced gap of YM-9 and the Casimir pinch of
YM-40. The strong-coupling Yang–Mills gap is the even part of the unit
turn, read on the lowest content.

## 4. What this says about the bridge

```text
                 odd part O                    even part E
order            first (the turn)              second (half the squared turn)
direction        circular (R)                  real (1)
protocol         deterministic refinement      symmetric independent records
law              weights: q, m                 squares: q², j(j+1)
physics          phase, quantum                heat, relaxation, mass gap at strong coupling
tie              E² − O² = I ,   memory = −O²
```

One arrow carries both. Observation — a symmetric record — deletes the odd
part and leaves its square as loss. Quantum and thermal are not two
theories joined by a constant; they are the two parities of the same
native arrow, tied by T1.

## 5. Certificate

Circular contents q = 1…6 at three rational turns: E, O by repeated
multiplication against the Chebyshev closed forms, T1, memory, the
character law, the bound by q², monotone approach along six shrinking
turns. su(2): T1 and memory on content ½ for one axis; record means on
degree one (4×4) and degree two (10×10, spectrum {1, χ_1/3} with
multiplicity nine); the three finite-turn rate formulas at r = 1, 1/4,
1/32. Six tests; a wrong character, a non-native turn and a one-sided
record are rejected.

## 6. Claim boundary

```text
E� − O² = I ;  RECORD MEMORY = −O²                                    PROVED (exact, finite turn)
CHARACTER LAW FOR CIRCULAR CONTENTS; RATE q²                          PROVED
CASIMIR LAW FOR su(2) CONTENTS ½, 1; 3/4 AND 5/4 RECOVERED            PROVED (free record model)
REAL UNIT = HALF THE SQUARED CIRCULAR UNIT, PER TURN                  PROVED as orders; NOT a physical calibration
k_BT IN TERMS OF ħ                                                    NO — needs one physical input: the turn rate
INTERACTING CASE (shared readings, weak coupling, continuum)          NOT TOUCHED
```

The relation between the two units is fixed per turn. To turn it into a
relation between physical constants one more number is needed: how many
turns per unit of physical time. That, and the interacting case, are the
open gates.

## 7. Reproduce

```text
python qc3_even_odd_arrow.py
python -m unittest test_qc3_exact
```
