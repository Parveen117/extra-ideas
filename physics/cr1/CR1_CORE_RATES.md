# CR1 — The pure numbers of the core: its two lowest rates in the sector unchanged by both turnings

Monty Dabas. 9 October 2026. Python 3.12, standard library only. Exact rationals; every count is an exact
elimination.

TC1 showed that the rates of the core of the constant modes are g^(2/3) times pure numbers. This stage computes
the first two of those numbers, for three and for four turns, and says exactly which of the statements about
them are certified.

Sources read before building: **TC1** (the three numbers, the layer, the scale), **TV1**; tools **DS1** (the
count of directions past a line), **AG1-A7** (a gap as a count); RKF **theorum/28** (an outward certificate
needs its own tail bound; without one a ladder certifies from one side only).

## The operator

C is a free matrix of three cut components and d directions, X = x₁x₂ + x₂x₃ + x₃x₁ in its squared singular
values, and

```text
h = −½·Δ + X ,        rates of −(g²/2)·Δ + X/g²  =  g^(2/3) × rates of h        (TC1-C6) .
```

Readings unchanged by turning the cuts and the directions are functions of the symmetric functions e₁, e₂, e₃
of the three numbers (e₁ = |C|², e₂ = X).

## Results

**R1 — the generator in three variables.** On such readings ½Δ is

```text
L = 3d·∂₁ + 2(d − 1)·e₁·∂₂ + (d − 2)·e₂·∂₃
    + 2e₁·∂₁₁ + 8e₂·∂₁₂ + 12e₃·∂₁₃ + (2e₁e₂ + 6e₃)·∂₂₂ + 8e₁e₃·∂₂₃ + 2e₂e₃·∂₃₃ .
```

The number of directions enters three coefficients and nothing else. Certified against the Laplacian in the
entries of C on e_i and e_i·e_j, which fixes a second-order operator: 9 polynomial identities each for d = 3
and 4, 8 for d = 5.

**R2 — free means by one rule.** Under the weight Exp(−w·e₁), for P of degree m (e₁ counts 1, e₂ 2, e₃ 3):

```text
mean(P) = mean(L·P) / (2·w·m) .
```

Every mean follows from mean(1) = 1. Equal to the means by pairing of the entries (d = 3, 4), and it returns
TC1's 6, 9, 3, 42, 72.

**R3 — the ladder, read by counting.** Readings m·Exp(−w·e₁/2), m a monomial of degree at most D; S and H their
exact matrices for 1 and h (67 × 67 at D = 10, symmetric exactly). The number of rates of the ladder below μ is
the number of negative directions of H − μ·S, counted by exact elimination; the k-th rate of h in this sector is
not above the k-th of any ladder.

```text
three turns      lowest rate < 5.1868        second rate < 8.0480
four turns       lowest rate < 8.0030        second rate < 11.5660

degree D     size     three turns: lowest, second       four turns: lowest, second
0            1        5.35782                           8.17778
2            4        5.20357    8.48939                8.01901    11.96276
4            11       5.18860    8.12261                8.00439    11.62342
6            23       5.18691    8.06006                8.00300    11.57419
8            41       5.18676    8.04857                8.00291    11.56663
10           67       5.1868     8.0467                 8.0029     11.5657
```

(Upper ends of exact brackets.) The ladder only falls as it grows. Its first rung is TC1's Gaussian reading.

**R4 — the lowest rate from below.** Two steps, both exact. (i) For one pair, as an operator in y with x held,
−(a/2)·Δ_y + b·|x × y|² ≥ |x|·√(2ab): the ground rate of TC1's layer. Every turn is the x of d − 1 pairs; giving
half of each turn's kinetic part to them leaves h ≥ Σ_turns [−¼·Δ + β·|c|] with β² = (d − 1)/2. (ii) For one
turn, the reading Exp(−k·|c|^(3/2)) has local rate (β − 9αk²/8)·r + (15αk/8)·r^(−1/2), whose least value by the
three-term mean is 3·(25αβ²/128)^(1/3). Hence

```text
(lowest rate)³  ≥  675·d³·(d − 1)/512 :        three turns  4.14 < lowest rate < 5.1868 ,
                                              four turns   6.32 < lowest rate < 8.0030 .
```

## The numbers

```text
                     lowest      second      difference (for the eye)
three turns          5.1867      8.047       2.86
four turns           8.0029      11.566      3.56
```

In the units of TC1 a rate of the constant modes is g^(2/3) times these, per size of the torus. The difference
is the distance between two upper bounds; it is not a certified gap.

## What is certified, and what is not

```text
the generator and the rule for means                         PROVED (polynomial identities)
upper bounds for the two lowest rates                        PROVED (exact counts on the ladder)
a window for the lowest rate                                 PROVED : [4.14, 5.1868] and [6.32, 8.0030]
a lower bound for the second rate, hence for the difference  NOT OBTAINED
```

What blocks the last line: step (i) of R4 throws away the quartic core and the correlation between pairs, and
by a rough estimate (not certified) the comparison it leaves has its second rate near 5.1 — under the first
rate of h itself. A bound for the
difference needs either an outward certificate for the ladder (theorum/28: a bound on what the readings above
degree D can still take away) or a comparison that keeps the core. Neither is built.

## What is not shown

- No mass gap: not for the core (previous lines), and not for the fabric.
- Only the sector unchanged by both turnings. Readings that turn under the directions (they are unchanged by
  the cuts, so they belong to the record) are not computed; one of them may lie below 8.05.
- The lowest rate is taken to exist as the bottom of a discrete set, and a positive reading with local rate at
  least λ is taken to bound it by λ (both named inputs, standard for this operator).
- In older terms (general knowledge, from memory, not re-read): with the usual normalization, which differs by
  2^(1/3), the two numbers for three turns are the known small-volume values 4.1167 and 6.386; the tests check
  that agreement to three places. The line's part: the generator in the three symmetric functions with d in
  three coefficients, the one rule for means, the reading of the ladder by counting, the elementary window for
  the lowest rate, and the statement of what a certified difference still needs.

## Claim boundary

```text
L IN e₁, e₂, e₃ ; d IN THREE COEFFICIENTS                                   PROVED (d = 3, 4, 5)
mean(P) = mean(LP)/(2wm)                                                    PROVED ; checked against pairing
LOWEST RATE < 5.1868, SECOND < 8.0480 (THREE) ; < 8.0030, < 11.5660 (FOUR)   PROVED (exact counts, degree 10)
(LOWEST RATE)³ ≥ 675 d³(d − 1)/512                                          PROVED (pair bound ; local rate)
DIFFERENCE 2.86, 3.56                                                       FOR THE EYE ; not a bound
SECOND RATE FROM BELOW ; OTHER SECTORS ; A GAP                              NOT OBTAINED ; NOT COMPUTED ; OPEN
```

## Reproduce

```text
python cr1_core_rates.py          (about half a minute)
python -m unittest test_cr1
```

## Later note (GC1, 9 October)

[GC1](../gc1/GC1_CORE_GAP_COUNT.md) counts with the comparison of R4 instead of taking only its lowest level:
under the row condition at most one reading lies under 5.5210 (three turns) and 8.0060 (four), so the second
rate has a floor and the gap is at least 0.3342 and 0.0031. With that, the lowest rate is in [5.1865, 5.1868]
and [7.9942, 8.0029] (the reading of this ladder at degree 10, its mean of h² exact). CM1's turning readings
close the gap from above: 2.3922 and 2.9796. Nothing above is changed.
