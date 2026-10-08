# DO1 — The observer is on the diagonal: seen = lost, and that alone gives the field of a centre

Monty Dabas. 8 October 2026. Python 3.12. Exact rational arithmetic and sympy.

The owner's correction: in the framework the observer is not fixed. The observer is on the diagonal, the line of
maximum information, and is inside the equation. SD1 and MM1 had fixed the cut along the block's own axes, where
nothing is lost, and stated the law as "no memory". This stage puts the observer where the owner says.

Sources read before building: RKF **theorum/24** Theorem 6.1 (S = R + D, F = R − D), **EMK-1** (determinant in a
plane), **DG1** (observed = lost at r = 2r_s; "a reason why the field is ρ = −½ from the diagonal alone: not
found"), **PH3-T2** (the native transport is the equal share of the two readings; share law 4λ(1−λ)), **QC4-T1**,
**LN1-N4/N5** (three cuts: the third layer returns the sign of b₁₂b₁₃b₂₃; the gravity element read across its
axis sees the same in both slots and has seen = lost), **SE1** (rates (ρ, 1, 1); e₂ = 0), **NC1** (ρ), **SD1**,
**MM1**, **GF1**, **CL1/HB1** (the eight-mark value 1/√2), **PT1-P3** (the quarter turn).

## 1. One plane

A block with rates λ₁, λ₂ read by cuts turned by χ from its own axes reads mean ± half-difference·cos 2χ and loses
(half-difference)²·sin² 2χ.

**P1.** On the diagonal, χ = 45°:

```text
both cuts read the same        a = mean
the lost part is largest       b² = (half-difference)²
F = R − D = a² − b²            is the determinant of the block (EMK-1), for each cut
```

Only there is the F of a single cut the determinant. And the one cut's pair (seen, lost) gives both rates,
λ = a ± b: the diagonal observer has the whole information. On the block's own axes nothing is lost and one cut
knows one rate.

So the determinant identity of EMK-1 is theorum/24's F = R − D for the observer on the diagonal. The observer's
own split is the equation.

## 2. Three cuts: the centre

Let the three cuts be equally inclined to the direction of fall. The stretch (ρ, 1, 1) of SE1 then reads

```text
every cut sees        c·(ρ + 2)/3
every pair loses      c·(ρ − 1)/3 .
```

**P2.** Seen = lost holds exactly for ρ = −½. In d cuts: ρ = −(d − 2)/2. Then every plane determinant vanishes
and each cut loses (d − 1) times what it sees. For three cuts the block is ½·[[1, −1, −1], [−1, 1, −1], [−1, −1, 1]].

This is SE1's law, NC1's ratio and MC1's exponent, from one statement: *the observer on the diagonal loses what
he sees.* DG1's open item is answered.

**P3 — the converse.** Take every block in which each cut reads a and each pair loses a, with any signs:

```text
two cuts       one kind:   rates (2, 0)·a                      a single reading; no field
three cuts     two kinds:  rates (3, 0, 0)·a                   a single reading; no field        sign of b₁₂b₁₃b₂₃ = +
                           rates (−1, 2, 2)·a                  the field of a centre             sign of b₁₂b₁₃b₂₃ = −
four cuts      three kinds, one with rates (1 − √5, 0, 2, 1 + √5)·a
```

In three cuts the diagonal observer with seen = lost sees either nothing or the field of a centre, and the one
sign that LN1-N4 says the third tower layer returns decides which. In four cuts a kind with irrational rates
appears.

## 3. Any stretch

**P4.** For any block there is a cut frame in which every cut reads the mean (built in two steps; the second is a
diagonal in one plane). In it

```text
e₂ = Σ_planes ( a² − K_ij² ) = Σ_planes F ,        empty space:   total lost = 2 × total seen .
```

SE1's "variance = (d − 1) × mean²" is "lost = (d − 1) × seen" for the diagonal observer. For a centre the balance
holds in every plane; in general it holds in total, and with a source Σ F is its count (CS1).

## 4. Records of turns

**P5.** For a record of turns the diagonal is the equal share of 1 and R, a quarter turn apart: seen = lost = ½,
and the mean has size 1/√2 — the eight-mark value of CL1 and HB1, and DG1's N = 1/√2 at r = 2r_s.

## What this settles

```text
the law, stated with the observer in it      F = R − D = 0 on the diagonal (seen = lost); with a source, F = its count
SD1 / MM1: "no memory (D = 0)"               the same invariant read on the block's own axes; the diagonal reading is
                                             the one in which the observer's split is the equation
DG1: "ρ = −½ from the diagonal alone"        FOUND: seen = lost for cuts equally inclined to the fall
three cuts                                   the only case with exactly: nothing, or a centre
```

## What is not shown

- P1 and P4 are the same invariants as before in another cut frame; no measured number changes.
- Why the observer is on the diagonal is the framework's statement (PH3: the native transport is the equal share);
  it is not derived here.
- P3 is an enumeration of sign patterns for two, three and four cuts; no meaning is claimed for the four-cut kinds.
- Beyond the stretch (the other components of the law, SW2's turning centre) the diagonal reading is not written.

## Claim boundary

```text
ON THE DIAGONAL OF A PLANE: EQUAL READINGS, LARGEST LOST PART, F = R − D = det                 PROVED (symbolic)
ONLY THERE IS A SINGLE CUT's F THE DETERMINANT ; ITS (SEEN, LOST) GIVES BOTH RATES             PROVED
CUTS EQUALLY INCLINED TO THE FALL: SEEN = LOST ⇔ ρ = −(d − 2)/2                                PROVED (d = 2 … 6)
THREE CUTS: EQUAL SEEN AND LOST ⇒ (3, 0, 0) OR (−1, 2, 2), BY THE SIGN OF b₁₂b₁₃b₂₃            PROVED (all sign patterns)
ANY BLOCK: e₂ = Σ_planes F FOR THE DIAGONAL OBSERVER ; LAW: LOST = 2 × SEEN                    PROVED
RECORD OF TURNS: DIAGONAL = QUARTER TURN, SIZE 1/√2                                           PROVED (exact)
WHY THE OBSERVER IS ON THE DIAGONAL                                                           THE FRAMEWORK'S STATEMENT
```

## Reproduce

```text
python do1_diagonal_observer.py
python -m unittest test_do1
```

## Later note (QD1, DU1, 9 October)

P1's turn to the diagonal is H = (C₁ + C₂)/√2, which exchanges the two cuts (DU1): for a two-valued record it is
the exchange of its two descriptions, and the record equal to its own diagonal reading has seen = lost. On the
diagonal the reading rule and the fixed-observer rule agree (QD1). Nothing above is changed.
