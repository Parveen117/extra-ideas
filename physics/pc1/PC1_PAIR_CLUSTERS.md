# PC1 — Pair clusters: a comparison that keeps the correlations; three turns from two, four from three, on all readings

Monty Dabas. 9 October 2026. Python 3.12, standard library only. Exact rationals; trial vectors are found in
floats and then used as exact rational vectors.

Every floor so far for the second rate of the core came from splitting it into single vectors — columns (GC1,
SG1) or rows (DR1) — and the split throws away the correlations between them: 15% of the lowest rate at three
turns. That is why the three-turn floor stood at 5.5210 against a lowest rate of 5.1867 and a ceiling of 7.58.
This stage splits the core into its **sub-cores** instead — three turns into its three pairs — and keeps, for
each sub-core, its own lowest reading and the floor behind it. The loss in the lowest rate is then 3%, and the
count needs no row condition.

```text
two turns        lowest rate 2.6592 … 2.6594        gap ≥ 0.5536      on all readings
three turns      second rate ≥ 5.6923 (was 5.5210)   gap ≥ 0.5056      on all readings (was 0.3342, physical only)
four turns       second rate ≥ 8.4699                gap ≥ 0.4670      on all readings (DR1: 0.4425, physical only)
```

Sources read before building: **CB1** (in full; merged here), **OM1** (the floor behind a reading; the overlap
bound (2)), **GC1** (G6: the core of n turns as the sum of its sub-cores; G4), **DR1** (floors ν(m) of one
vector), **CR1**, **CM1** (§5: a count needs an actual floor on everything left out); RKF **theorum/28**.

## The step

GC1-G6: the core of n turns is the sum of its n sub-cores of n − 1 turns, each in the form a·T + b·X with
a = 1/(n − 1), b = 1/(n − 2), and a dilation makes each one κ = (a²b)^(1/3) times the core h of n − 1 turns.
G6 used only the lowest rate E of h. Suppose more is known: on **all** readings of n − 1 turns, h has its
lowest reading g and nothing else under a floor ρ. Then h ≥ ρ − (ρ − E)·P_g, and

```text
P1      h_n   ≥   κ·[ n·ρ − (ρ − E)·M ] ,        M = Σ over the n sub-cores of  P_g ⊗ 1 .
```

**P2 — the levels of M.** M = J·J*, where J sends n one-turn readings (a₁ … a_n) to Σ_k g ⊗ a_k (g on the
turns other than k). The other product is small and explicit:

```text
J*·J  =  1 + R ⊗ (ones − 1) ,
```

R the one-turn part of g — a positive operator of trace 1, with levels r₁ ≥ r₂ ≥ … summing to 1 — and "ones"
the n × n table of ones. So the levels of M are 1 + (n − 1)·r on equal parts and 1 − r on parts that sum to
zero. Only one of them is above 1 + (n − 1)(1 − r₁).

**P3 — the count.** Hence at most **one** reading of h_n, among all readings, lies under

```text
Z_n  =  κ·[ (n − 1)·ρ + E − (n − 1)(1 − r₁)(ρ − E) ] ,
```

and r₁ ≥ ⟨a ⊗ … ⊗ a, g⟩² for any one-turn reading a. With r₁ = 1 (a lowest reading without correlations)
this is κ[(n − 1)ρ + E]: one sub-core raised to its floor. The correlations enter only through 1 − r₁, and
r₁ is large: 0.956 for two turns.

The step needs a floor on all readings of the sub-core, with no row condition: a sub-core is cut out of the
core, and what is unchanged by the row turning of all n turns is not unchanged by that of n − 1.

## Two turns: the start

By columns, with halves (GC1): h₂ ≥ k(c₁) + k(c₂), unit 1/2. Sort all readings by their row turning and by
their parity in each turn; the comparison keeps all three.

```text
no row turning, even in both turns        one reading at the bottom ; the next at (e₀ + e₁)/2 = 3.2130
row turning two, even in both             3.2921
row turning four ; two turns with turning number two            4.0652 ; 4.2461
odd in a turn                             3.3592
```

The last line is the new one. With halves, a reading odd in one turn has floor 2.8487 — the cheap reading
that only the row condition removed in GC1. But for a sector's lowest level no count is needed, so the split
need not be even: give **all** of the other turn's kinetic part to the layer. Then
h₂ ≥ −½·Δ₁ + √2·|c₁|, whose unit is 1, and on readings odd in c₁ it is at least ν(5) = 3.3592 (DR1).

So on all readings of two turns at most one lies under ρ₂ = 3.2130. With one exact reading of CR1's ladder
(degree 14 in the two numbers e₁, e₂; the third is zero for two turns) and its mean of h²:

```text
lowest rate of two turns  in  [ 2.6592 , 2.6594 ] ,        gap ≥ 0.5536 on all readings .
```

Its one-turn share: the product reading a ⊗ a with a = (1 + (2/25)|c|²)·Exp(−w|c|²/2) has squared overlap
0.9596 with the ladder reading (exact; its average over the directions is 1 + b·e₁ + b²(e₁²/8 + e₂/2)), the
ladder reading has squared overlap at least 0.9998 with the true lowest reading (OM1's bound), and angles
add: r₁ ≥ 0.9556.

## Three turns, and four

```text
three turns      κ = 2^(−2/3)        Z₃ = κ·[ 2·3.2130 + 2.6592 − 2·0.0444·(3.2130 − 2.6592) ]  =  5.6923
```

At most one reading of the three-turn core lies under 5.6923, among all readings. With GC1's reading
(5.1868) the gap is at least 0.5056; with CM1's turning reading, at most 2.3922. The lowest rate, by GC1-G4
with the new line: [5.1865, 5.1868].

The same step again, now with the three-turn core as the sub-core (floor 5.6923 on all readings, lowest rate
from both sides, one-turn share at least 0.9564 from the Gaussian product):

```text
four turns       κ = 18^(−1/3)       Z₄ = 8.4699 ,        gap ≥ 0.4670 on all readings .
```

Two controls. Without the full layer for odd turns the two-turn floor is 2.8487 and Z₃ falls back under
5.5210. And knowing nothing of the one-turn share (r₁ = 1/2) does the same: both inputs are needed.

## What holds the number down

Z₃ is κ·(2ρ₂ + E₂) less a small term, so every 0.1 added to the two-turn floor ρ₂ adds 0.126 to the
three-turn line. The floor 3.2130 is the comparison's radial step of one turn, in the sector with no row
turning and both turns even. In floats (a rough basis, not certified) the second level of two turns on all readings appears
near 3.66, a six-fold level of the kind odd in a turn, with the even sector's own second level above it; with
3.66 the same step would give 6.23. The comparison by columns cannot say this: its second level in that sector is 3.2130
whatever is done. What is needed is a floor for the second even level of **two** turns — a problem in three
variables, smaller than the one this stage started from.

## What the certificate checks

```text
P1    the shares of turns and pairs ; κ³ = a²b
P2    J*J = 1 + R ⊗ (ones − 1) on exact small models (2 and 3 points, 3 and 4 turns) ; trace R = 1
two   the sector floors and their least ; the two units ; the control with halves ; the reading and its
      window ; the gap ; the share r₁
three Z₃ against 5.5210 ; the gap ; two controls ; the lowest rate with the new line
four  Z₄ and the gap
```

Tests: the cluster bound on small models where it is attained (three turns of three points, plane rotations
in floats); the direction average of the product reading against exact pairing of the entries; monotone
behaviour of the share and the line.

## Inputs taken from outside

```text
as in GC1 and DR1: one vector in its radius and its sphere ; levels of a self-adjoint operator and counts
a positive operator of trace 1 has levels summing to 1 ; J J* and J* J have the same levels away from zero
the angle between two readings is at most the sum of their angles to a third
CR1's ladder ; DR1's floors ν(5), ν(7), ν(11) ; CM1's turning reading
```

## What this changes in the line

```text
GC1, SG1, DR1    second rate of three turns       was: ≥ 5.5210 (physical readings)      now: ≥ 5.6923 (all readings)
                 gap of three turns               was: ≥ 0.3342                          now: ≥ 0.5056
DR1              gap of four turns                was: ≥ 0.4425 (physical)               now: ≥ 0.4670 (all readings)
GC1-G6           two turns                        was: lowest rate ≤ 2.7445              now: 2.6592 … 2.6594, gap ≥ 0.5536
CM1 §5           a floor on everything left out   the cluster step gives one on all readings, turn after turn
```

The earlier notes are not rewritten; later notes in GC1 and DR1 point here.

## What is not shown

- **The floor is still far from the ceiling**: 0.51 against 2.39 for three turns. The whole distance is now in
  one number, the second even level of two turns.
- **No volume.** The same step can be written for a lattice — the operator as a sum of its plaquette pieces,
  each with its lowest reading and floor. But the count compares a line with the lowest rate itself, and the
  distance between the two grows with the number of pieces while the gap does not; I expect it to close only
  for a few pieces, as it does here for a few turns (not tried). CB1's bound for a lattice still falls with
  volume, and nothing here changes that.
- Which reading is the first excited one on all readings of three or four turns; ceilings on all readings.
- No mass gap; no fabric.
- In older terms: sums of sub-system operators for lowest rates are a known device, and the levels of a sum
  of projectors through its small product are standard; I did not search for this count in print. The line's
  part: the floor on all readings of two turns by the uneven layer, the count through the one-turn share, and
  the numbers.

## Claim boundary

```text
h_n ≥ κ[nρ − (ρ − E)M] ; LEVELS OF M ARE 1 + (n − 1)r AND 1 − r                        PROVED
AT MOST ONE READING OF TWO TURNS UNDER 3.2130, ON ALL READINGS                         PROVED, given the inputs
TWO TURNS: LOWEST RATE IN [2.6592, 2.6594] ; GAP ≥ 0.5536 ; r₁ ≥ 0.9556                 PROVED, given the inputs
THREE TURNS: AT MOST ONE READING UNDER 5.6923, ON ALL READINGS ; GAP ≥ 0.5056           PROVED, given the inputs
FOUR TURNS: AT MOST ONE READING UNDER 8.4699, ON ALL READINGS ; GAP ≥ 0.4670            PROVED, given the inputs
THE SECOND EVEN LEVEL OF TWO TURNS (FLOATS: THE STEP WOULD THEN GIVE 6.23)              NOT OBTAINED
ANY STATEMENT WITH VOLUME ; A MASS GAP                                                 OPEN
```

## Reproduce

```text
python pc1_pair_clusters.py          (about three seconds)
python -m unittest test_pc1
```
