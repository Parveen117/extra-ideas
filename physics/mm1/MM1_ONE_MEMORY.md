# MM1 — One quadratic form for every memory: e₂. Gravity's law is "no memory between the cuts", and TP1's coefficients follow from two no-memory conditions

Monty Dabas. 8 October 2026. Python 3.12. Exact rational arithmetic and sympy, with the functions of CV1 and TP1.

OR1 left two premises: the record rule, and equivalence of local frames; and it found gravity's field law to be
a ratio (SE1), not a least variance. The question put to this stage: is that ratio the record rule among the three
cuts?

Sources read before building: **EMK-1** (the determinant of a block in a plane: seen product − lost square),
**GE2-T2** (record mean S, memory I − S†S; T5: a deterministic record has none), **IN1-T1/T2/T4/T5** (reading
tensor ρ = ½(n + r·C); det ρ = ¼(n² − r·r); a single reading is null; the invariant is frame-independent; a record
of two readings has n² − r·r = 2p(1−p)(n₁n₂ − r₁·r₂), the unrecoverable memory), **QC4-T1** (a single reading has
no variance), **LN1-N1/N2** ((H²)_ii = a_i² + lost_i; every tower level carries it), **TP1-T2/T4** (equivalence
selects 1 : 2 : −4; the profile r_s/r is kept on a plane of coefficients), **SE1**, **OR1**, **EG1** (source: the
energy reading; constant k not derived), **SW1-T2** (a constant swirl is no field), **GW1** (a local turn of the
frame is not a wave).

## 1. The form

e₂(X) = ½[(tr X)² − tr X²]: for a block in a plane it is the determinant (EMK-1); for three cuts it is the sum of
the three plane determinants.

**U5 — what it measures.** In any cuts, with a_i = X_ii seen and lost_i = (X²)_ii − a_i² (LN1):

```text
2·e₂(X)  =  ( Σ a_i )²  −  Σ ( a_i² + lost_i )  ,        tr T_λ(X) − T_λ(tr X) = −2λ·e₂(X) .
```

e₂ is the failure of the uncut reading (the trace) to commute with the first generation of the tower.

**U1 — records.** For records X_s with weights of sum one,

```text
mean of e₂  −  e₂ of the mean   =   Σ_{s<t} p_s p_t · e₂( X_s − X_t ) .
```

## 2. Three memories, one form

```text
turns (GE2-T2, OR1)         e₂(u) = 1 for a turn ;  memory = 1 − e₂(S)                         U2
readings (IN1)              e₂(ρ) = 0 for a single reading ;  e₂(record) = ¼(n² − r·r)        U3
                            = ½ p(1−p)(n₁n₂ − r₁·r₂) : the unrecoverable memory, the rest mass²
stretch of the frame (SE1)  e₂(P) = 0 for a stretch along a single cut ;                       U4
                            e₂(K) = Σ_{i<j} λ_i λ_j (1 − overlap²) for K = Σ λ_i P_i
```

The record memory of turns, the unrecoverable memory of readings and the quantity in gravity's law are the same
quadratic form. The sign of e₂ on a difference is + for turns and − for self-dagger readings: a record of turns
loses e₂, a record of readings gains it. In both, (uncut)² = (held by the cuts) + 2e₂.

## 3. Gravity's law as a no-memory law

By SE1, on flat slices G₀₀ = e₂(K). So:

```text
empty space      e₂(K) = 0      the stretch, read through the cuts, has no unrecoverable part:
                                it behaves as a single reading
with a source    e₂(K) = k·u    the unrecoverable memory of the stretch is the energy of the source (EG1)
```

**U6 — the coefficients.** On flat slices TP1's three invariants give

```text
Q = A·tr K² + B·(turn)² + Γ·(tr K)² ,        A = 2a₁ + a₂ ,   B = 2a₁ − a₂ ,   Γ = a₃ .
```

Ask of Q what holds for every memory above:

```text
a pure turn of the frame has no memory        (GE2-T5: a deterministic record; SW1-T2)       B = 0
a stretch along a single cut has no memory    (IN1-T2, QC4-T1: a single reading has none)    A + Γ = 0
```

These two conditions leave a₁ : a₂ : a₃ = 1 : 2 : −4, TP1's law, and nothing else. The second condition alone is
exactly TP1-T4's plane, 2a₁ + a₂ + a₃ = 0: the profile r_s/r is kept precisely when a single cut has no memory.

**U7 — in time.** For every fall velocity on flat slices the volume V of a cell of the falling frame obeys

```text
V″ / V  =  2·e₂(K) − R₀₀ ,
```

so in empty space it is linear in the frame's own time. For the centre at rest the cell labelled r₀ has volume
factor (r₀^{3/2} − (3/2)·√r_s·t)·√r₀, exactly.

## What this settles

```text
OR1: "gravity's field law is a different kind"      on flat slices it is the same kind: no memory (e₂) between the cuts,
                                                     equal to the source's energy where there is one
TP1: coefficients from equivalence                   on flat slices also from two no-memory conditions
TP1-T4: the plane that keeps r_s/r                   = "a single cut has no memory"
the premises                                         on flat slices one: no memory in a shared record
```

## What is put in

- SE1's family (one time, flat slices) and TP1's three invariants as the available quadratic laws.
- The reading of a stretch along one cut as a single reading, and of a turn of the frame as a deterministic
  record. These are the identifications; the algebra then fixes the coefficients.

## What is not shown

- Beyond flat slices. In general frames TP1-T2's equivalence argument remains the proof of 1 : 2 : −4; that the
  two no-memory conditions suffice there is not shown.
- U1–U5 are elementary algebra of the second symmetric function. The line's part is the identification of three
  of its certified quantities (GE2's, IN1's, SE1's) as that one form, and U6.
- Why the unrecoverable memory of the stretch equals the *energy* of the source, and the constant k, are not
  derived (EG1's open item).
- Only the energy component of the law is read this way; the other components are not.
- The record rule itself remains a premise.

## Claim boundary

```text
2e₂ = (Σ SEEN)² − Σ(SEEN² + LOST) ; TRACE AND TOWER COMMUTE UP TO e₂                PROVED
MEAN OF e₂ − e₂ OF THE MEAN = Σ p_s p_t e₂(X_s − X_t)                              PROVED (symbolic)
GE2's TURN MEMORY, IN1's INVARIANT AND SE1's LAW ARE THE FORM e₂                    PROVED (exact instances; identities)
ON FLAT SLICES: NO MEMORY FOR A PURE TURN AND FOR A SINGLE CUT ⇒ 1 : 2 : −4        PROVED
SINGLE-CUT CONDITION = TP1-T4's PLANE                                              PROVED (against TP1's stored result)
V″/V = 2e₂ − R₀₀ ; LINEAR VOLUME IN EMPTY SPACE                                    PROVED (symbolic)
THE SAME FOR GENERAL FRAMES                                                        NOT SHOWN
WHY MEMORY OF THE STRETCH = ENERGY OF THE SOURCE ; THE CONSTANT                    NOT DERIVED
THE RECORD RULE                                                                    PREMISE
```

## Reproduce

```text
python mm1_one_memory.py
python -m unittest test_mm1
```
