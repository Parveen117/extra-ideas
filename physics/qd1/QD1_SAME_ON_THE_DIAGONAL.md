# QD1 — On the diagonal the reading rule and the fixed-observer rule are the same

Monty Dabas. 9 October 2026. Python 3.12. Exact counting, sympy, mpmath.

The owner's statement: on the diagonal quantum and classical are the same. This stage makes it exact in three
records of the line and finds where the two part.

Sources read before building: **T24-6.1** (S = R + D, F = R − D), **IN1** (a reading read by a cut; lost to one
cut = observed by the others), **SD1** (4R·D = the other cuts' F²), **DO1** (the observer is on the diagonal),
**DG1** (sin 45° = cos 45°; share law q(1 − q), largest ¼ at ½), **CL1/HB1** (the eight-mark value 1/√2),
**RD1** (thermal light: mean count x/(1 − x)), **PT1-P3** (the quarter turn).

"Classical" is read here as the rule with the observer fixed outside: a mark set beforehand, read by its sign.
"Quantum" is the line's reading rule: seen cos²(θ/2), lost sin²(θ/2).

## Results

**D1 — one reading, one cut.** Seen + lost = 1, flip F = cos θ, spread = 4·seen·lost. On the diagonal the cut
reads a fair coin with the largest spread: exactly what it reads from a record with no direction at all. One
cut on the diagonal cannot tell the single reading from the fully mixed record; the difference is with the
other cuts (IN1).

**D2 — a pair, two rules.** For two cuts at angle θ read on a pair:

```text
reading rule            flip = cos θ
fixed-observer rule     flip = 1 − 2θ/π        (exact count on clocks of 8, 16, 32, 64 marks: 1 − 4k/Q)
```

They agree on the cut (θ = 0, π) and on the diagonal (θ = π/2), and nowhere else. Any rule that changes sign
when seen and lost are exchanged is zero on the diagonal: there every rule gives the same reading. That is the
owner's statement, exact.

On the four-mark clock (cuts and diagonals only) every sum over pairs is the same for both rules. They first
part on the eight-mark clock: 1/√2 against ½ — the eight-mark value of CL1 and DG1. For the sum over four pairs:

```text
fixed observer (every assignment, by enumeration)      2
reading rule                                           2√2 = 2.82843
measured (arXiv:1506.01865)                            2.82759 ± 0.00051
```

The measurement is 1.6 errors below the reading rule and 1600 errors above the fixed observer. The owner's
earlier correction — the observer is not fixed outside — is what this number decides.

**D3 — counts.** For a mode with mean count n, the spread of the count is n(1 + ηn):

```text
repeating record (p_k ∝ x^k)        n + n²        S = R + D
marks placed independently          n             R
two-valued record (empty or full)   n − n²        F = R − D          with R = n , D = n²
```

The three kinds of counting are the three members of T24-6.1 for one split. The two-valued spread is the share
law q(1 − q), largest ¼ at ½.

**D4 — the diagonal mode.** R = D at n = 1, x = ½: the mode is empty half the time. Its quantum is kT·ln 2:

```text
energy  1 bit × kT ,     entropy  2 bits ,     free energy  −1 bit × kT .
```

The curvature of the entropy, −1/s″, is n for the counting law, n² for the wave law and n + n² for the actual
law: on the diagonal the two classical pictures weigh the same.

**D5 — numbers.** For the sky's thermal light (2.7255 K) the diagonal mode is at 39.4 GHz. Below it — on the
wave side — lie 7.9 % of the count and 1.3 % of the energy.

## What is put in

- A pair whose members read opposite on every common cut; the line has the single reading (IN1), not this pair.
- The reading of "classical" as "observer fixed outside". The temperature of the sky is recalled.

## What is not shown

- D2's two rules and the sum 2√2 are known (general knowledge); D3 is the known spread of the three statistics.
  The line's part: they agree exactly on the diagonal and for a stated reason; the three spreads are S, R, F of one
  split; the place where the rules part is the eight-mark value.
- Why the split of D3 is R = n, D = n² is not derived from T24's lift; it has the form.

## Claim boundary

```text
RULES AGREE ON THE CUT AND THE DIAGONAL, NOWHERE ELSE ; EXCHANGE-ODD RULES VANISH ON THE DIAGONAL    PROVED
FIXED OBSERVER: SUM ≤ 2 (ENUMERATION) ; READING RULE: 2√2 ; EIGHT-MARK VALUE 1/√2 AGAINST ½          PROVED ; KNOWN NUMBERS
MEASURED 2.82759 ± 0.00051                                                                          MATCH (reading rule)
SPREADS n + n², n, n − n² = S, R, F ; DIAGONAL MODE AT kT·ln 2                                       PROVED (form)
```

## Reproduce

```text
python qd1_same_on_the_diagonal.py
python -m unittest test_qd1
```
