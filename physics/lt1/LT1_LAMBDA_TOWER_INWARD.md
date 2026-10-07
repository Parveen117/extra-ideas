# LT1 — The λ-tower on the gravity field: inward, between two references, to a single reading

Monty Dabas. 8 October 2026. Python 3.12, standard library only. Exact
rational arithmetic.

The owner's statement: raising λ in the tower moves inward in the
diagram, where more is held; any λ can be taken as the reference; and
at infinite λ the secret of counting may be hidden.

Sources read before building: **RMG2** (the λ-tower
T_λ(H) = H + λH²: radial on the disc; ℓ′ = ℓ + log((1+λλ₊)/(1+λλ₋)),
between ℓ and 2ℓ; fixed sets — the centre and the nilpotent boundary;
T1: the number sectors), **NC1** (response element of the gravity
field, tanh²η = m = r_s/r), **GR1** (N = sech η), **IN1** (the
invariant), **SY1** (dual sector), **MO1**.

## Results

**T1 — one generation moves inward.** At infinite λ a generation is
squaring: the rapidity doubles, m′ = 4m/(1+m)², and the place moves

```text
r′ = ( r + r_s )² / 4r .
```

Inward for every r > r_s; fixed at the horizon and at infinity. Far
away r′ ≈ r/4. Iterating, every place goes to the horizon:
40 → 10.5 → 3.15 → 1.37 → 1.02 → … (in units of r_s).

**T2 — λ moves continuously.** For finite λ the eigenvalue ratio of the
element goes from q² (λ = 0) to q⁴ (λ = ∞), increasing all the way:
q³(1+λq)/(q+λ). Changing λ is a continuous move inward.

**T3 — two references, one generation apart.** With r = 4s the map of
T1 is

```text
r′ = s ( 1 + r_s/4s )² ,
```

the relation between the radius s in which space is conformally flat
and the radius r′ of MO1 in which space is flat. In each of them the
law has the same form — memory = r_s/radius — at its own level of the
tower:

```text
level below      tanh²(η/2) = r_s / 4s        N = (1 − x)/(1 + x),  x = r_s/4s
level of MO1     tanh² η    = r_s / r′
```

So "stopping at a level and calling it flat" is exact here: the two
standard ways of laying out space around a mass are two consecutive
levels of the λ-tower.

**T4 — what the tower removes.** The invariant of the normalised
element is 1 − ρ². A generation sends it to ((1−ρ²)/(1+ρ²))², and
iterating sends it to zero. At the end the element is rank one: a
single null reading (IN1). The horizon is where the response element of
the field has become one reading.

**T5 — at the end, powers count.** At that boundary the quarter-turned
element squares to zero (RMG2-T1), and there (1 + tN)^k = 1 + k·tN
exactly: composition is addition, and the count is exact.

## Reading

The owner's three sentences, as certified statements:

```text
raising λ moves inward                         T1, T2
any level can be the reference                 T3 — and two levels are the two known layouts of space
the end of the tower is where counting lives   T4, T5 — the element becomes a single reading and powers add
```

## What is not shown

- T5 is a fact about the algebra at the boundary. That the counting of
  WQ1 (hν per wave) *comes from* this boundary is not shown; no link
  between the two was derived.
- Whether a count belongs to the horizon itself — a whole number of
  units of something there — is a known conjecture in physics and is
  not derived here.
- T3 identifies two levels with two layouts of space; levels above and
  below these two have not been given a meaning.
- "More is held inward" is, in this stage, the memory m rising to one
  and the invariant falling to zero; no density was defined.

## Claim boundary

```text
GENERATION: m′ = 4m/(1+m)², r′ = (r+r_s)²/4r, INWARD, HORIZON FIXED            PROVED
FINITE λ INTERPOLATES MONOTONICALLY                                           PROVED
ONE GENERATION = CONFORMALLY FLAT RADIUS → FLAT-SPACE RADIUS                  PROVED
INVARIANT → 0; BOUNDARY ELEMENT RANK ONE                                      PROVED
EXACT COUNTING IN THE BOUNDARY SECTOR                                         PROVED (SY1)
COUNTING OF WAVES OR OF THE HORIZON DERIVED FROM THE BOUNDARY                 NOT SHOWN
```

## Reproduce

```text
python lt1_lambda_tower_inward.py
python -m unittest test_lt1_exact
```
