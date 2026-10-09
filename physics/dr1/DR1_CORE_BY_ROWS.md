# DR1 — The core read by rows: three cuts in d directions; a gap for three to ten turns, and the end of the walk in the number of turns

Monty Dabas. 9 October 2026. Python 3.12, standard library only. Exact rationals; levels and trial vectors are
located in floats and then certified exactly.

GC1 and SG1 compared the core by its columns: d turns, each with three cut components. That comparison gives
0.3342 for three turns, 0.0031 for four, and nothing from five turns on. OM1 then showed how a floor behind a
computed reading returns to it. This stage reads the same weight by its rows — three cuts, each with d
components — and the count changes: every d from three to ten gets a gap, four turns moves from 0.0031 to
0.4425, and the walk in the number of turns of GC1-G6 gets its end point in closed form.

Sources read before building: **OM1** (in full; merged here), **SG1**, **GC1** (G1 the sign certificates, G2
the count, G4 the two-sided lowest rate, G6 the walk in d), **CR1** (the ladder for any d), **TC1** (C1: the
weight is one number of the three squared singular values, the same by rows and by columns), **CM1**; RKF
**theorum/28** (a floor of its own before any count).

## The weight read two ways

C has three rows r₁, r₂, r₃ (each a point of d directions) and d columns c₁ … c_d (each a point of three
cut components). The weight of the core is

```text
X  =  Σ_{i<j} |c_i × c_j|²  =  Σ_{a<b} ( |r_a|²·|r_b|² − (r_a·r_b)² ) ,
```

the second symmetric function of the three squared singular values (TC1), which do not care whether the
matrix is read by rows or by columns. The row turning — the gauge condition — mixes the three rows; turning
the d directions turns every row in its own space.

## Results

**D1 — the comparison by rows.** Give half of each row's kinetic part to its two pairs. For a pair (a, b) the
row r_b is an oscillator in the d − 1 directions across r_a, of frequency |r_a|/2, so its ground rate is
(d − 1)·|r_a|/4. What is left:

```text
h   ≥   k′₁ + k′₂ + k′₃ ,        k′ = −¼·Δ_d + ((d − 1)/2)·|r| ,
```

three copies of one operator in d directions. By columns it was d copies in three directions with
β = √((d − 1)/2); at d = 3 the two are the same.

**D2 — one row.** On readings of turning number l in d directions, with |r| = ℓ·s,

```text
k′  =  unit × ( −u″ + [(m − 1)(m − 3)/(4s²)]·u + s·u ) ,       m = d + 2l ,      unit³ = (d − 1)²/16 :
```

a turning number l in d directions is the level problem of d + 2l directions. Let ν(m) be its first level and
ν′(d) the second level for l = 0. The solution is u = s^p·w(s), p = (m − 1)/2, with w an exact power series;
where the potential is still above λ it is positive and rising, and from there on the equation bounds it step
by step (grid 1/64), as in GC1-G1. Certified floors (each 0.002 to 0.004 under the level):

```text
m        3        4        5        6        7        8        9        10       11       12       13       14
ν(m)     2.3381   2.8700   3.3592   3.8155   4.2461   4.6561   5.0489   5.4270   5.7924   6.1466   6.4910   6.8264
ν′(m)    4.0879   4.4910   4.8824   5.2609   5.6277   5.9837   6.3301   6.6676
```

In closed form, from the reading s^p·Exp(−k·s^(3/2)) with k² = 2/9:   ν(m) ≥ (3/4)·(2m − 1)^(2/3).

**D3 — the count.** The row turning contains the three half turns that flip two rows at once. A reading
unchanged by them has its three rows all even or all odd. The comparison keeps the parity of each row, and
its lowest reading of one row is even. So among physical readings:

```text
all rows even          one reading at the bottom ; any other has a row above the bottom, at least
                       2ν(d) + min( ν′(d), ν(d + 4) )
all rows odd           at least 3ν(d + 2)
```

and at most **one** physical reading of h lies under

```text
z′ = unit × min( 2ν(d) + min(ν′(d), ν(d + 4)) , 3ν(d + 2) ) .
```

Only the three half turns are used, not the whole row turning. By columns the cheap reading was one turn with
turning number one, removed by the whole row turning; by rows it is one row with turning number one, removed
by parity.

**D4 — the gap for three to ten turns.** From above, one exact rational reading of CR1's ladder (degree 10
for three and four turns, 6 for the others), with its mean of h² exact; from below, the lowest rate by GC1-G4
with the line z′:

```text
d      line z′     lowest rate              gap ≥      ρ_d³ = E³/(d³(d − 1))
3      5.5210      5.1865 … 5.1868          0.3342     2.58363 … 2.58399
4      8.4454      8.0028 … 8.0029          0.4425     2.66950 … 2.66957
5      11.6008     11.0664 … 11.0707        0.5301     2.71052 … 2.71361
6      14.9597     14.3543 … 14.3573        0.6024     2.73856 … 2.74028
7      18.5022     17.8361 … 17.8384        0.6639     2.75711 … 2.75815
8      22.2126     21.4927 … 21.4948        0.7179     2.77019 … 2.77095
9      26.0776     25.3096 … 25.3115        0.7662     2.77998 … 2.78056
10     30.0858     29.2745 … 29.2762        0.8097     2.78759 … 2.78804
```

For three turns rows and columns agree (the same 5.5210). For four turns the gap floor is 0.4425 where the
columns gave 0.0031, and with CM1's turning reading the gap is in [0.4425, 2.9710]. From five turns on the
columns give nothing — their line is under the lowest rate — and the rows give the first floors. The windows of
ρ_d³ rise one after another, as GC1-G6 says they must.

**D5 — the end of the walk.** Three copies of the closed floor of D2 from below, the Gaussian reading from
above:

```text
(729/1024)·(d − 1)²·(2d − 1)²   ≤   E(d)³   ≤   (729/256)·d³·(d − 1)        for every d ≥ 3 ,
```

and the left side is at least (1 − 2/d) times the right. So ρ_d = E(d)/(d·(d − 1)^(1/3)), which never falls
(GC1-G6), rises to exactly

```text
ρ_∞  =  (729/256)^(1/3)  =  9/(2·32^(1/3))  =  1.41741… ,          0 ≤ 1 − ρ_d³/ρ_∞³ ≤ 2/d .
```

With many turns the lowest rate of the core is the Gaussian one: the correlations between turns that cost
3.2% at three turns and 2.1% at four die out like 1/d. (At d = 3 the left side is CR1's 4.14.)

## What the certificate checks

```text
D1    the weight by rows and by columns on exact matrices (d = 3, 4, 5, 7) ; the oscillator of a pair
D2    the kept polynomials against the equation ; tails ; 11 first-level and 7 second-level sign certificates ;
      four refusals 0.02 above ; the order of the floors ; the closed floor
D4    eight exact readings (⟨ψ, hψ⟩ two ways) ; line above reading for d = 3 … 10 ; the lowest rate from both
      sides ; rows = columns at three turns ; four turns against 0.0031 ; five to ten turns against the
      columns ; the gap window of four turns ; the order of the ρ windows
D5    the closed bounds, their ratio (d ≤ 400), and the certified windows between them
```

## Inputs taken from outside

```text
the kinetic part of one row in its radius and its sphere of d − 1 dimensions, with levels l(l + d − 2) on
    the sphere and parity (−1)^l
levels of a self-adjoint operator and the count under one condition (as in GC1)
CM1's turning reading for four turns ; CR1's ladder
```

## What this changes in the line

```text
GC1 / SG1    gap of four turns                 was: 0.0031 / 0.0030           now: 0.4425
GC1          five or more turns                was: nothing                    now: floors for d = 5 … 10
GC1-G6       the limit of ρ_d                  was: not shown, ≤ 1.4175        now: exactly (729/256)^(1/3)
GC1-G4, OM1  lowest rate of four turns         [7.9942, 8.0029], [8.002732, 8.002896]     now: [8.00283, 8.00290]
OM1          hidden floor behind the four-turn reading (δ = z − v/(z − η))     8.0052 with z = 8.0060 ; with z′ = 8.4454 it is 8.4453
```

The earlier notes are not rewritten; a later note in each points here.

## What is not shown

- **Three turns — the case of three directions of space — gains nothing**: there the rows are the columns, and
  the floor stays 0.3342 against a ceiling of 2.39. Both comparisons leave out the correlations between
  the vectors they split, and the cheapest excitation they allow is one vector's radial step.
- **Every d.** Beyond ten turns nothing is certified. In floats the margin keeps growing, like 0.37·d^(1/3);
  the two branches that have closed floors (a row with turning number two; all rows odd) stay above the
  Gaussian value for every d tried up to 20000, but the radial second level ν′(d) has no closed floor here.
- No turning reading from above for five or more turns, so no ceiling for their gaps.
- No mass gap, no fabric, no volume: as in GC1.
- In older terms: vectors-with-many-components limits of this kind are usually had from a mean-field
  argument; I did not search whether the closed bounds of D5 or the count by rows are in print. The line's
  part: the same weight read by rows, the parity count from three half turns, exact floors for one row in
  any number of directions, and the end point of GC1's walk.

## Claim boundary

```text
X BY ROWS = X BY COLUMNS ; h ≥ k′₁ + k′₂ + k′₃                                           PROVED
FLOORS OF ν(m), m ≤ 14, AND ν′(d), d ≤ 10 ; ν(m) ≥ (3/4)(2m − 1)^(2/3)                   PROVED (exact series ; local rate)
AT MOST ONE PHYSICAL READING UNDER z′ ; GAP > 0 FOR d = 3 … 10 (TABLE)                   PROVED, given the inputs
LOWEST RATES OF d = 3 … 10 FROM BOTH SIDES (TABLE)                                      PROVED, given the inputs
(729/1024)(d − 1)²(2d − 1)² ≤ E(d)³ ≤ (729/256)d³(d − 1) ; ρ_d → (729/256)^(1/3)          PROVED, given the inputs
A GAP FOR EVERY d ; A SHARPER FLOOR AT THREE TURNS                                      NOT OBTAINED
A MASS GAP ; THE FABRIC ; ANY VOLUME                                                    OPEN
```

## Reproduce

```text
python dr1_core_by_rows.py          (about seven seconds)
python -m unittest test_dr1
```
