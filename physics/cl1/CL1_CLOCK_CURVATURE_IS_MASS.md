# CL1 — The native clock answers the clock question: mass² is clock curvature

Monty Dabas. 7 October 2026. Python 3.12, exact rational arithmetic only.

PR4 ended by asking what fixes a position-dependent unit of time N and
could not derive a law for it. Owner's direction: time is already derived
in this repository — clock-free calculus, recognition-seam calculus —
look at the mathematics first and interpret afterwards.

Sources read before building:

- **T24** (Recognition-Kernel-Framework theorum/24): arrows, cuts and
  composition residues with no time parameter; "a real parameter may
  later provide one chart of the arrow system; it is not the clock of the
  calculus".
- **GE2-T3**: the history ledger τ = Σ ½|v|² — additive on histories, not a
  state function, not a unique clock.
- **R39** (echo clock), **R41** (relational clock): **R41.1** process
  relative to the tick record; **R41.2** residue of a closed history,
  Q⁻¹‖(W^Q − I)w‖²; **R41.3** tick–phase curvature D·T·D†·T† = ζ_Q,
  [D,T]†[D,T] = 2 − ζ_Q − ζ̄_Q, and [N, P] = ιR; **R41.4** the
  count–response inequality.
- **R38**: c_Σ = 1/√2; one Z block contains eight source events.
- **PR1** (coin-turn walk), **QC3-T1** (E² − O² = I).

## 1. What the mathematics says

**T1 (clock curvature).** For a tick T and a phase readout D with tick
phase ζ = c + ιs: D·T·D⁻¹·T⁻¹ = ζ and the clock curvature is 2 − 2c.
(R41.3, here for every rational turn, on marks without a wrap.)

**T2 (mass² is clock curvature).** The coin-turn walk of PR1 with the
same ζ has, on its uniform sector, 2 − (U + U⁻¹) = 2 − 2c. So

```text
mass²  =  clock curvature  =  2 − ζ − ζ̄ ,          cone speed  =  1 − (clock curvature)/2 .
```

A sector with rest turn θ is a clock with tick phase θ; what PR1 called
its mass is the curvature R41 derived for the tick against its phase.

**The dyadic ladder** (R38's positive half-roots):

```text
marks Q      tick phase     clock curvature     cone speed
   2            π               4                 −1         R41: record phase and successor are H and K
   4            π/2             2                  0         nothing propagates
   8            π/4             2 − √2             1/√2      R28's coin; R38's c_Σ; eight source events per block
   Q → ∞        → 0             → 0                → 1       light-like
```

R38's sharp speed is the speed of the eight-mark clock.

**T3 (count against the tick).** With E_T = (T + T⁻¹)/2, O_T = (T − T⁻¹)/2
and the count N (away from a wrap):

```text
[N, O_T] = E_T ,          [N, E_T] = O_T .
```

The count turns the tick's odd part into its even part and back. R41.3's
[N, P] = ιR is this, with P = ι·O_T and the wrap written out.

**T4 (a coin that varies from place to place).** For coins (c_x, s_x):

```text
(U + U⁻¹)ψ₁(x) = c_{x−1} ψ₁(x−1) + c_x ψ₁(x+1) + (s_x − s_{x−1}) ψ₂(x−1) ,
(U + U⁻¹)ψ₂(x) = c_{x+1} ψ₂(x+1) + c_x ψ₂(x−1) + (s_{x+1} − s_x) ψ₁(x+1) .
```

The speed becomes local, and the two components are exchanged by exactly
the difference of the odd part between neighbours. A constant coin has
no exchange term. This is the native form of PR3-T4's term gρ·K.

**T5 (history residue).** For Q ticks of a turn W = Exp(φR):
Q⁻¹‖(W^Q − I)w‖² = 2(1 − E(Qφ))/Q. The history closes exactly when
W^Q = I: the quarter turn closes at Q = 4, 8; the rational turn
(3 + 4R)/5 never does (checked to Q = 8). This is R41.2 for a turn, and
it is QC1's silence condition qΘ ∈ 2πℤ.

## 2. The answer to PR4's question

- **There is no position-dependent unit of time in the native
  calculus.** Time is a count of events along a history (T24, GE2-T3,
  R41.1). The function N of PR3–PR4 is a chart, not a native object, and
  "what fixes N″" is not a native question. PR4's refusal stands, for a
  sharper reason.
- **The native question it turns into has an answer.** The curvature a
  clock carries is the curvature of its tick against its phase, and that
  is the sector's mass² (T2). What can vary from place to place is the
  coin, and when it does the exact new term is the difference of the odd
  part between neighbours (T4).
- **Rates differ between histories, not between places.** Two histories
  with the same ends carry different counts and different ledgers; the
  failure of a closed history to return is the residue of T5.

## 3. Physics reading (after the mathematics)

- A massive sector is a clock, and its mass is how curved that clock is.
  A light-like sector is a clock with no curvature: it cannot tick
  against its own phase.
- The speed of a sector and the curvature of its clock are one number:
  speed = 1 − curvature/2. A sector propagates as fast as its clock is
  flat.
- The factor ρ of PR3 is then a comparison of counts between two
  histories in relative boost, not a field on space.
- A coin that changes with position acts as an exchange between the two
  light-like readings proportional to the local change of mass.

## 4. Certificate

T1–T2 for five rational tick phases (including 1 and ι): Weyl relation on
six marks, curvature, rest defect of the walk on both components, the
speed law. The ladder at Q = 2, 4 exactly and Q = 8 through its square.
T3 as matrix commutators on a window of seven marks. T4 on a line of
seven sites with seven different rational coins, against the constant
coin. T5 for two turns, Q = 1…8. Seven tests; a non-turn and a corrupted
walk are rejected.

## 5. Claim boundary

```text
MASS² = CLOCK CURVATURE = 2 − ζ − ζ̄ ; SPEED = 1 − CURVATURE/2             PROVED
DYADIC LADDER; R38's SPEED AS THE EIGHT-MARK CLOCK                         PROVED for Q = 2, 4, 8 (the last through its square)
[N, O_T] = E_T , [N, E_T] = O_T                                           PROVED away from a wrap
VARYING COIN: EXCHANGE = DIFFERENCE OF THE ODD PART                        PROVED
HISTORY RESIDUE; CLOSURE ⇔ W^Q = I                                         PROVED for turns
A FIELD EQUATION FOR A UNIT OF TIME                                        REFUSED — not a native object
WHICH Q NATURE USES; ANY PHYSICAL MASS OR SPEED                            NO
GRAVITY                                                                    NOT TOUCHED
```

R41's own results are cited, not re-derived; T1 and T3 restate R41.3 in
the even/odd language on marks without a wrap.

## 6. Reproduce

```text
python cl1_clock_curvature_is_mass.py
python -m unittest test_cl1_exact
```
