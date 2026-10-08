# FD1 — One function for every memory: det(1 + zX)

Monty Dabas. 9 October 2026. Python 3.12. sympy and exact rational arithmetic.

MM1: every memory of the line is the quadratic form e₂. SH1: the N-fold memory of a record is e_N. The owner's
question: can the determinant of the whole record, in our style, carry the line? This stage builds
det(1 + zX) = Σ zⁿ eₙ(X) on the cut carrier and reads the line's results off it.

Sources read before building: **T24-6.1** (S = R + D, F = R − D, Gram forms), **EMK-1** (determinant of a block),
**IN1** (reading tensor; det ρ = ¼(n² − r·r)), **MM1** (e₂), **SH1-S4** (e_N = Πp × Gram determinant; exclusion is
independence), **SD1** (seen × lost = the other cuts' F²), **SE1, DO1** (the field of a centre: ρ = −(d − 2)/2),
**OR1, RC1** (closing circuits), **QD1-D3** (three kinds of count), **DU1** (the two-valued chain), **theorum/28**
§4, §9 (declared tails, outward certificate), **YM-6** (inertia across a cut).

## Results

**F1 — the function.** det(1 + zX) = 1 + z·e₁ + z²·e₂ + … ; for IN1's reading tensor it is
1 + z·n + z²·¼(n² − r·r). The reciprocal function 1/det(1 − zX) carries the partner numbers hₙ, and at second
order

```text
h₂ = R + D ,     e₂ = R − D ,     R = ½ (tr X)² ,   D = ½ tr X² ,     exp(z·tr X) has R alone.
```

That is T24-6.1 with S = h₂ and F = e₂: the determinant is the F side, its reciprocal the S side. The vacuum
law of gravity (SE1, DO1) is "the z² coefficient is zero": for the stretch (ρ, 1, …, 1) it gives ρ = −(d − 2)/2.

**F2 — closed circuits.** −log det(1 − zK) = Σ zⁿ tr Kⁿ / n: the function is built from closed circuits only,
weighted by length. For a record of marks and allowed steps, tr Kⁿ counts closed n-step walks and

```text
1 / det(1 − zK)  =  Π over primitive circuits  1 / (1 − z^length)
```

(checked by enumeration to length 9). Closure, the line's rule since RC1 and OR1, is what the function is made of.

**F3 — across a cut.** With P + Q = 1,

```text
det(1 − zK) = det(1 − z·seen block) × det(1 − z·lost block − z²·exchange·(1 − z·seen block)⁻¹·exchange) .
```

For a Gram form at second order: e₂ = e₂(seen) + e₂(lost) + seen·lost − exchange². For a single reading e₂ = 0,
so seen × lost = exchange²: SD1's identity, now one coefficient of one function.

**F4 — counts through a cut.** Take N independent readings as one exclusive record (SH1). The number of them a
cut P sees has the law

```text
Σ s^k P(k seen)  =  det( 1 + (s − 1)·K_P ) ,      K_P the seen kernel,
```

so the count is that of independent coins whose biases are the readings of K_P: mean tr K_P, spread
tr K_P(1 − K_P). **Seen through any cut, the exclusive record counts like classical coins; on the diagonal they
are fair coins.** The three kinds of count of QD1 are one form, det(1 − ηwN)^(−η): the determinant, its
reciprocal, and between them exp(w·tr N).

**F5 — rates.** The zeros of det(1 − zB) are the inverse readings; a rate is the log-ratio of two zeros. For DU1's
chain it is 2k*. The diagonal turn leaves the function unchanged, and doubling the cell is
det(1 − z²B²) = det(1 − zB)·det(1 + zB).

**F6 — an infinite ladder.** For readings qⁱ, i = 0, 1, 2, …: eₙ = q^(n(n−1)/2) / ((1 − q)…(1 − qⁿ)), and the
ladder cut at M steps misses at most the fraction q^(M−n+1)/(1 − q) of each coefficient. A declared tail in the
shape of theorum/28 §4; no limit is taken without it.

## What this gives

The line's separate identities — IN1's invariant, MM1's memory, SH1's exclusion, SD1's seen × lost, the three
counts, the rate of a chain — are coefficients, factors or zeros of one function. ONE_LAW's "one form e₂"
becomes "one function, read at second order".

## What is not shown

- Every identity here is a known one for finite matrices (general knowledge). The line's part is the placement.
- F4 puts in SH1's rule that an exclusive record is the squared wedge of its readings.
- F6 is one ladder. For a general infinite record the tail has to be declared and bounded case by case; nothing
  is claimed without it.
- A gap in an infinite record is a statement about where the zeros are. The function restates that question;
  it does not answer it. The weak-coupling wall of MG1 stands.

## Claim boundary

```text
det(1 + zX) = Σ zⁿeₙ ; h₂ = R + D, e₂ = R − D ; VACUUM LAW = SECOND COEFFICIENT ZERO           PROVED
−log det = SUM OVER CLOSED CIRCUITS ; PRODUCT OVER PRIMITIVE CIRCUITS                          PROVED (enumeration to length 9)
SPLIT ACROSS A CUT ; SINGLE READING: SEEN × LOST = EXCHANGE²                                   PROVED
COUNT THROUGH A CUT = det(1 + (s − 1)K_P) ; FAIR COINS ON THE DIAGONAL                         PROVED (exact, three sizes)
RATE = LOG-RATIO OF ZEROS ; LADDER qⁱ WITH DECLARED TAIL                                       PROVED
LOCATION OF ZEROS FOR THE COLOUR RECORD                                                        OPEN
```

## Reproduce

```text
python fd1_one_function_for_every_memory.py
python -m unittest test_fd1
```
