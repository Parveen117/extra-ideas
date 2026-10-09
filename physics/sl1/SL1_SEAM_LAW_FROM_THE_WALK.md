# SL1 — The seam law from the walk: the mirror t → 1/t is proved from the diagonal step, two means, and one count at the weak end

Monty Dabas. 9 October 2026. Python 3.12, standard library only. Exact integer series (order 400 in q), lattice
enumerations, directed rational enclosures (tools/pm1, AG1).

The largest borrowed statement of this line is the seam law of the pair record — S(1/t) = t·S(t) and its two
companions. AG1 held it as "agreement certified; equality classical, pinned", and TW1, HL1 and PM1 lean on it.
The owner's rule is that proofs use the line's own steps: move along the diagonal from state to state, without
an outside reference. This stage proves the seam law that way.

Sources read before building: **AG1** (A1 the seam law as pinned; A3, A4 the diagonal step and the mean),
**HL1** (H2 sheets; H3 the clock M(S, L) halves), **TW1**; tools **PM1** (P5), **DS1**; Publications **LAM-1**
(the theta seam as enclosures).

## The record

q = Exp(−πu), and for whole n

```text
a = Σ q^(n²) ,   b = Σ (−1)ⁿ q^(n²) ,   c = Σ q^((n+½)²) ;        S = a² ,  F = b² ,  L = c² .
```

A prime is one diagonal step, q → q² (u → 2u). M(x, y) is the common limit of the repeated step
(x, y) → ((x + y)/2, √(xy)); it lies between the two numbers at every step and M(λx, λy) = λ·M(x, y).

## The proof

**M1 — one step, three lattice identities.**

```text
a² = a′² + c′² ,        b² = a′² − c′² ,        c² = 2·a′·c′ .
```

Each is one substitution. Pairs (m, n) with m + n even are (x + y, x − y), and m² + n² = 2(x² + y²). Pairs with
m + n odd are (j + k + 1, j − k), and m² + n² = 2(j + ½)² + 2(k + ½)². On the half lattice
(m + ½)² + (n + ½)² = (x² + y²)/2 with x + y odd. So

```text
S = S′ + L′ ,      F = S′ − L′ ,      L² = 4·S′·L′ ,      and from these  S² = F² + L² .
```

(Series to order 400; the substitutions on 625 and 400 points.) The triangle of AG1-A2 is a consequence here,
not an input.

**M2 — two means.** One step of the mean on (S, F) is the pair at 2u — (S + F)/2 = S′ and S·F = F′² — so
M(S, F) = lim S(2ⁿu) = 1 (AG1-A4). One step of the mean on (S, L) is, read from M1 at u/2,

```text
( (S + L)/2 , √(S·L) ) = ½·( S(u/2) , L(u/2) ) ,
```

half the pair at half the heat time. Repeating, and because the mean lies between its two numbers,

```text
2⁻ⁿ·L(u/2ⁿ)  ≤  M(S, L)(u)  ≤  2⁻ⁿ·S(u/2ⁿ)        for every n .
```

At u = 1 the two ends are 0.8346, 1.1803 at n = 0; 0.99254, 1.00748 at n = 1; equal to ten digits at n = 3.

**M3 — the count at the weak end.** √t·a(t) is the sum of Exp(−πx²) over the points x = n·√t with step √t, and
√t·c(t) the same over the half-shifted points. The weight falls on each side of zero, so each sum is caught
between two areas:

```text
1 − √t ≤ √t·a(t) ≤ 1 + √t ,        1 − √t ≤ √t·c(t) ≤ 1 + 2√t        (the area under Exp(−πx²) is 1) .
```

Hence (1 − √t)² ≤ t·L(t) and t·S(t) ≤ (1 + √t)². Put t = u/2ⁿ in M2 and multiply by u:

```text
(1 − √t)²  ≤  u·M(S, L)(u)  ≤  (1 + √t)²        for every t = u/2ⁿ ,        so        u·M(S, L)(u) = 1 .
```

The record reads its own heat time: 1/u is the mean of its source and its lost part. (HL1-H3 had the halving
and the value as agreement; here it is proved.)

**M4 — a right triangle is fixed by its two means.** Let x² = y² + z² with M(x, y) = 1 and M(x, z) = w. Write
y = x·k′, z = x·k. Then M(1, k)/M(1, k′) = w, and the left side rises strictly with k (M(1, ·) rises); so k,
and then x = 1/M(1, k′), are fixed. Now take w = u:

```text
( S, F, L ) at 1/u         M(S, F) = 1 (M2) ,   M(S, L) = u (M3 at 1/u)
( uS, uL, uF ) at u        M(uS, uL) = u·M(S, L) = 1 (M3) ,   M(uS, uF) = u·M(S, F) = u (M2)
```

Both are right triangles with means 1 and u. They are the same triangle:

```text
S(1/u) = u·S(u) ,        F(1/u) = u·L(u) ,        L(1/u) = u·F(u) ,
```

and for one turn, by positive roots, a(1/u) = √u·a(u) and b(1/u) = √u·c(u). ∎

The mirror exchanges the flip and the lost part because it exchanges the two means.

## The one input

The area under Exp(−πx²) is 1. This is where π enters and it is the only thing taken from outside the walk.
Control: with the weight Exp(−3u) in place of Exp(−πu) every step above is the same, M(S, F) is still 1, and
u·M(S, L) comes out π/3 = 1.0472 (enclosure, 100 digits) — the count at the weak end is the one place that
knows the constant.

## What the certificate checks

```text
M1    the three identities and their consequences as series ; the substitutions on boxes
M2    the bracket 2⁻ⁿL(u/2ⁿ) ≤ M(S, L) ≤ 2⁻ⁿS(u/2ⁿ) at u = 1, n = 0 … 6, closing in ; M(S, F) = 1
M3    the two inequalities at t = 1/4, 1/16, 1/64, 1/256 ; u·M(S, L) = 1 at u = 1/4, 1, 3 (100 digits)
M4    M(1, k) rising ; the seam law at u = 2, 3/2, 5, 1/3 with both triples' means ; one turn at u = 4, 9/4, 1/9
```

The inequalities of M3 are far from sharp: √t·a(t) − 1 is 7·10⁻⁶ at t = 1/4, 3·10⁻²² at 1/16, 10⁻⁸⁷ at 1/64.
The proof needs only that they close.

## What this changes in the line

```text
AG1-A1    seam law of the pair                              was: agreement, classical        now: PROVED (one input)
AG1-A2    S² = F² + L²                                      was: series                      now: from M1's substitutions
HL1-H3    u·M(S, L) = 1                                     was: agreement                   now: PROVED
AG1-A6, TW1-W4, HL1-H5, HL1-H6                              mirror description of the turn block, the limit π/2,
                                                            the first sheet term, the seam values: each is one
                                                            change in u of the law above (term by term)
PM1-P5    seam law of the prime-turn series                 k = 0 is this law ; k ≥ 1 still agreement
```

The earlier notes are not rewritten; a later note in each points here.

## What is not shown

- The area under Exp(−πx²) is taken, not derived.
- PM1's series with the weights α^(4k), k ≥ 1, and anything about the Riemann line beyond its theta seam.
- Nothing about a mass gap.
- In older terms this is the inversion law of the theta functions, usually had from a summation formula; that
  it follows from Landen's step and the arithmetic–geometric mean is, I expect, known (not re-read). The line's
  part: the proof by its own steps — the diagonal substitution, the two means, the squeeze from the weak end,
  the triangle fixed by its means — with the single outside input named and tested by a control.

## Claim boundary

```text
a² = a′² + c′² ; b² = a′² − c′² ; c² = 2a′c′ ; S² = F² + L²                  PROVED (substitution ; series)
2⁻ⁿL(u/2ⁿ) ≤ M(S, L)(u) ≤ 2⁻ⁿS(u/2ⁿ) ; M(S, F) = 1                           PROVED
u·M(S, L)(u) = 1                                                            PROVED, given the area under Exp(−πx²)
S(1/u) = uS ; F(1/u) = uL ; L(1/u) = uF ; a(1/u) = √u·a ; b(1/u) = √u·c       PROVED, given the same
WITH Exp(−cu): u·M(S, L) = π/c                                               CERTIFIED at c = 3 (control)
THE AREA UNDER Exp(−πx²) ; THE SERIES WITH α^(4k), k ≥ 1                      TAKEN ; NOT SHOWN
```

## Reproduce

```text
python sl1_seam_law_from_the_walk.py
python -m unittest test_sl1
```

## Later note (SG1, 9 October)

[SG1](../sg1/SG1_GAUGE_SINGLET_GAP.md), section 6, points out that the line "Hence (1 − √t)² ≤ t·L(t)" in M3
squares 1 − √t ≤ √t·c(t), which is allowed only when 1 − √t ≥ 0, that is for t ≤ 1; for larger t the lower end
is zero. The proof uses the line at t = u/2ⁿ as n grows, where t ≤ 1, so M3 and everything after it stand as
written. The correction is to the range of that one line. Nothing else above is changed.
