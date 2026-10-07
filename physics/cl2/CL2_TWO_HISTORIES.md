# CL2 — Two histories with the same ends: the turned one counts less, and the cost is a split-sector defect

Monty Dabas. 7 October 2026. Python 3.12, exact rational arithmetic only.
Continues CL1.

CL1 concluded that in the native calculus rates differ between
histories, not between places. This stage computes the difference
exactly: two histories between the same two events, one straight, one
turned.

Sources read before building: **R41.1** (a process advances by the count
of its tick), **GE2-T3** ("a path and its reversed point motion can close
at the endpoint while retaining a positive history ledger"), **PR2-T1/T2**
(the form X² = x² + y² − t², the boosted mass generator), **CL1-T2**
(the rest turn is the clock), **QC3-T3** (the character law).

## 1. The count of a leg

A leg is a displacement X = xH + yK + tR with X² = −τ² on the upper
sheet. A sector of rest turn g carried along it has energy–momentum
P = (g/τ)X (PR2-T2).

**T1.** The count of the leg is −sc(P·X) = g·τ: the rest turn times the
form-length of the displacement.

## 2. The corner cost

**T2.** For two legs X₁, X₂ and the straight history X₁ + X₂ of
form-length T,

```text
T² − (τ₁ + τ₂)² = 2 τ₁ τ₂ ( cosh δ − 1 ) ,        cosh δ = −B(X₁, X₂) / (τ₁ τ₂) ≥ 1 ,
```

with B the bilinear form. Equality holds exactly when the legs are
collinear. The straight history has the largest count; every turn of the
velocity costs count, and the cost is the defect (E − 1) of the boost
between the legs — the split-sector twin of QC3's record defect (1 − E).

Example: legs (x, t) = (3, 5) and (−3, 5), each of form-length 4. Straight
history: 10. Turned history: 8. Corner: cosh δ = 17/8.

## 3. Many corners: the split character

**T3.** N equal legs on a hyperbola of radius ρ, corner boost 2b:

```text
count / straight count = N / χ_N(b) ,          χ_N(b) = sinh(N b) / sinh(b) ≥ N ,
```

with equality only for N = 1. This is the character of QC3-T3 in the
split sector:

```text
circular sector (turns, records)     χ_q(θ) = sin(qθ)/sin(θ) ≤ q      the record loses phase
split sector   (boosts, histories)   χ_N(b) = sinh(Nb)/sinh(b) ≥ N    the turned history loses count
```

For fixed ends, a finer turned history has a smaller count
(total rapidity 8a with cosh a = 41/40):

```text
legs      1        2          4         8
ratio     1        0.7015     0.6370    0.6215        (exact fractions in the result file)
```

The smooth limit — the uniformly accelerated history of PR3, count g·ρ·η —
is the infimum of this sequence.

## 4. Both sheets

**T4.** Under g → −g every count changes sign and no ratio changes. Mass
and antimass carried along the same two histories disagree by the same
factor.

## 5. Reading

- The native form of "acceleration slows a clock" is a statement about
  counts on two histories; no unit-of-time field appears.
- What is lost is fixed by the corners alone: 2τ₁τ₂(cosh δ − 1) per
  corner, in squares. A history with no corner loses nothing.
- GE2-T3 already says that a history can close at its endpoint and still
  retain a positive ledger. T2 is that statement with its exact value in
  the boost sector.
- The same character governs what a record loses (QC3) and what a turned
  history loses (here); the sector decides the direction of the
  inequality.

## 6. Certificate

T1 on nine rational legs (two space directions). T2 on all 36 pairs —
cost formula, cosh δ ≥ 1, equality exactly on the one collinear pair —
and a three-leg history. T3 on a hyperbola with cosh a = 41/40: points by
the addition law, every chord, the four polygons, the identity
count·χ_N = N·straight, χ_N > N, and the monotone order. T4 by the sheet
flip. Six tests; a spacelike leg and a non-boost are rejected.

## 7. Claim boundary

```text
COUNT OF A LEG = g τ                                                     PROVED
CORNER COST 2τ₁τ₂(cosh δ − 1); STRAIGHT HISTORY HAS THE LARGEST COUNT    PROVED
POLYGON RATIO N/χ_N; SPLIT CHARACTER ≥ N; FINER ⇒ SMALLER                PROVED on the certified family
SMOOTH LIMIT η / sinh η                                                  STATED (infimum; not certified as a limit)
SAME RATIOS ON BOTH SHEETS                                               PROVED
THE WALK (LATTICE) VERSION: A MOVING PACKET'S OWN PHASE                  NOT BUILT
GRAVITY: HISTORIES IN A REGION WHERE THE COIN VARIES (CL1-T4)            NOT TOUCHED
```

## 8. Reproduce

```text
python cl2_two_histories.py
python -m unittest test_cl2_exact
```
