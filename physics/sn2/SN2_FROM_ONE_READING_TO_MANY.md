# SN2 — The spread of a count from one reading to many: 1/3 → 1/4; and what adds along one chain

Monty Dabas. 9 October 2026. Python 3.12. sympy, exact enumeration over products of cuts, mpmath.

SN1 gave 1/4, 1/2, 1/3 with the even spread of the sector's angle put in. This stage asks what the line itself
can say about that spread.

Sources read before building: **SN1**, **FD1-F4** (count through a cut; seen kernel), **IN1** (Σ F_i² = 1: lost to
one cut = observed by the others; three cuts complete), **DO1** (no cut preferred), **EMK-1** (cuts C₁ = K,
C₂ = RK, C₃ = ιR), **DU1-U3** (joining cells), **DG1-D4** (memory = tanh²).

## Results

**N1 — what the number is.** spread / mean = Σq(1 − q)/Σq: the mean of the lost share, weighted by the seen
share. SN1's three sectors, read so:

```text
boost    the seen-weighted amplitude tanh η is evenly spread on 0 … 1 :   mean of amplitude² = 1/3
shear    the turn angle is evenly spread :                                mean of sin² = 1/2
turn     the same with the seen weight cos² :                             1/4
```

**N2 — one reading, no cut preferred, d cuts.** Since Σ F_i² = 1 (IN1), ⟨F²⟩ = 1/d and

```text
spread / mean = (d − 1) / 2d :      d = 2 : 1/4 ,    d = 3 : 1/3 ,    d → ∞ : 1/2 .
```

The flip F of such a reading is evenly spread only in three cuts, and ⟨F²⟩ equals spread/mean only in three cuts.
A single cut-complex reading (three cuts, IN1) with no preferred direction has 1/3.

**N3 — many readings.** Take n factors (carrier 2ⁿ), the cut P = ½(1 + C₁ on the first factor), and for the
readings' frame every Π = ½(1 ± σ), σ any product of cuts other than 1, none preferred. By enumeration
(n = 1, 2, 3) and by a closed count for every n, with N = 2^(n−1) readings:

```text
mean = N/2 ,      spread / mean = N² / (4N² − 1) :      1/3 ,  4/15 ,  16/63 ,  …  →  1/4 .
```

One reading gives the three-cut number 1/3; many readings give the turn sector's 1/4. The two ends of SN1's
table are one formula. No continuous average is used: the frames are the finitely many products of cuts.

**N4 — one chain: what adds.** Two cells with boost angles η₁, η₂ joined through a turn φ:

```text
1/seen = cosh²η₁ · cosh²η₂ · |1 + tanh η₁ · tanh η₂ · e^(iφ)|² .
```

With no turn preferred, the mean of log(seen) is exactly the sum of the two: **log(seen) adds**. The angle does
not: the mean of 1/seen is cosh²η₁cosh²η₂(1 + tanh²η₁tanh²η₂), not cosh²(η₁ + η₂).

## Verdict on SN1's put-in

```text
turn sector's even angle (many readings)       N3 gives its number 1/4 as the limit of an exact finite count
three-cut 1/3                                  N2, N3 at one reading
boost sector's even angle from joining cells   REFUSED for one chain: log(seen) adds, the seen share of a long
                                               chain falls off, and the spread is not even. The even spread is a
                                               property of many readings together; not derived here.
```

## What is not shown

- N² /(4N² − 1) and (d − 1)/2d are known numbers for these averages (general knowledge). That the average over
  products of cuts equals the average over all frames at this order is known and is not proved here; N3 states
  the finite average itself.
- N4's adding of log(seen) is the known behaviour of one disordered chain.
- The next number of N2's model is 0 in every d; SN1's sectors give 0, 1/4, 1/15. The two routes to 1/3 agree at
  this order only.

## Claim boundary

```text
SPREAD/MEAN = SEEN-WEIGHTED MEAN OF LOST ; ONE READING IN d CUTS: (d − 1)/2d                          PROVED
N = 2^(n−1) READINGS OVER PRODUCTS OF CUTS: N²/(4N² − 1) ; 1/3 → 1/4                                 PROVED (enumeration n ≤ 3; count for all n)
ONE CHAIN: MEAN LOG(SEEN) ADDS ; THE ANGLE DOES NOT                                                  PROVED
EVEN SPREAD OF THE BOOST ANGLE                                                                       NOT DERIVED
```

## Reproduce

```text
python sn2_from_one_reading_to_many.py
python -m unittest test_sn2
```
