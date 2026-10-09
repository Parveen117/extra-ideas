# MG1 — A gap as a mass count: the place of the proton's number, and why its value is still open

Monty Dabas. 9 October 2026. Python 3.12. Symbolic algebra (sympy), numbers with mpmath.

AC1's later note: one rule is missing — what fixes the mass count n = r_s/2ℓ of a kind of matter. The framework's
Yang–Mills line certifies gaps. This stage writes a gap as a mass count and reads off what the certified gaps give.

Sources read before building: **AC1** (later note), **BR1** (rest rate g, record area κ), **SW2/GM1** (lengths of a
centre). Yang–Mills line: **YM-1** (reduced gap −Log(I₂(2)/I₁(2))), **YM-9, YM-22** (free rate 3/4 at every
cutoff), **YM-21** (chain: decay exactly geometric, rate −Log r, r = I₂(κ)/I₁(κ), every length), **YM-40** (content
ladder), **YM93, YM97** (gap windows for two and three colours), **YM98** (block reduction at every coupling; §5
the weak-coupling obstruction), **YM99** (no uniform ladder in the quadratic sector).

## Results

**M1 — three lengths, one area.** For a centre with the turn of one reading, mass length × turn length is the
record area: r_s·a = ℓ². So

```text
mass count  n = r_s / 2ℓ = ℓ / 2a = (rest rate) × ℓ .
```

**M2 — a gap is a mass count.** The gap γ of a generator, per cell of length a_c, is the rest rate of its lightest
record: n = (ℓ / a_c)·γ. This is the place of a mass count in the Yang–Mills line.

**M3 — what is certified.** In cell units:

```text
YM-1   one graph, κ = 2                        0.8367
YM-9, YM-22   free rate, two colours            3/4      at every cutoff
YM95, YM98    free rate, three colours          4/3
YM93, YM97    windows                           ≥ 16/125 ,  ≥ 2/5
YM98          window ends                       > 0.37 ,  > 0.78
```

Every one is of order one. So the lightest colour record has a mass count of about ℓ/a_c: with the cell at the
record length, a mass near the unit itself. The proton's count is 7.69 × 10⁻²⁰. An order-one rate reaches it only
with a cell of 2.8 × 10⁻¹⁶ m — which is putting the proton's size in.

**M4 — the one rate known at every coupling.** On the chain the rate is −Log(I₂(κ)/I₁(κ)), decreasing, and at
large κ it is 3/(2κ): a power. The proton's count would need κ = 2 × 10¹⁹. A large number in for a large number
out: refused as an explanation.

**M5 — what kind of law would do.** log(1/n) is 44.0 for the proton and 51.5 for the electron. A rate that falls
as an exponential of the coupling needs a number near 44, not 10²⁰.

## Verdict

```text
place of a mass count in the Yang–Mills line        n = (ℓ / a_c) × gap                    FOUND
value of the proton's count                         NOT GIVEN by any certified statement:
                                                    certified rates are of order one; the exact chain rate is a power
where it would have to come from                    the gap at weak coupling in four dimensions — the open wall of
                                                    YM98 §5 and YM99 (the compact, non-quadratic part of the vacuum)
```

The physics line's missing rule and the Yang–Mills line's open wall are the same place.

## Known physics, outside the line

In the standard theory the smallness of the proton's count is read as an exponential of the inverse coupling
(leading order, three colours, no matter: exponent 1/(2b₀g²), b₀ = 11/16π²). Taken at the record length this
gives g² ≈ 0.16, about θ ≈ 2 × 10² in the Yang–Mills line's coupling, against certified windows θ ≤ 1/64. This is
general knowledge with the prefactor and higher orders dropped; it is not a result of the line and nothing above
rests on it. It says only how far the wall is: a number of order 10², not 10²⁰.

## What is put in

- The reading "gap = rest rate of the lightest record" (M2).
- The constants of M5 are recalled, not re-read at source.

## Claim boundary

```text
r_s·a = ℓ² ; n = r_s/2ℓ = ℓ/2a = REST RATE × ℓ                                       PROVED
CERTIFIED YANG–MILLS RATES ARE OF ORDER ONE IN CELL UNITS                            READ FROM YM-1 … YM98
CHAIN RATE → 3/(2κ)                                                                  PROVED (numerically to κ = 10⁷)
THE PROTON'S MASS COUNT                                                              NOT DERIVED ; place found
```

## Reproduce

```text
python mg1_mass_count_from_a_gap.py
python -m unittest test_mg1
```

## Later note (DU1, 9 October)

M5 asked for a law that falls as an exponential of the coupling. DU1 finds one natively for a chain of two-valued
marks: rate = 2k* ≈ 2e^(−2k), with k → k − ½ ln 2 on doubling the cell. The chain of turns of M4 stays a power.
The proton's count is still not derived. Nothing above is changed.

## Later note (CT1, 9 October)

M4's κ = 2 × 10¹⁹ is k_c = 22.4 in the coupling of the chain's centre record (CT1): the same fact. Which of the two
is the number nature sets is not decided by the line. Nothing above is changed.

## Later note (CZ1, 9 October)

On a closed surface the centre record closes at fixed cosets, and the weight of a twist is an exponential of the
turn coupling: for a cube 4.27·κ·e^(−(6 − 3√3)κ), against the power 3/(4κ) for one open face. A twist weight of
7.7 × 10⁻²⁰ is κ = 61.7. This is not a mass count. Nothing above is changed.

## Later note (TW1, 9 October)

On the one non-abelian record with exact data (a closed surface, heat weight) the weak end is reached by walking
back from the strong side with a loss bounded once; no expansion at the weak end is used (TW1-W4, W5). The wall
of the verdict — the same for the four-dimensional fabric — stands. No mass count is derived. Nothing above is
changed.

## Later note (TV1, 9 October)

For the constant modes of d turns the quadratic layer leaves the weight (Σ sin²α_i)^−(d−2) on the valley:
closed for three turns, a logarithm at the centre points for four (TV1). This places where the non-quadratic
part of the verdict's wall must enter; it does not pass the wall. Nothing above is changed.
