# PT1 — The primes inside the native algebra: exact turns, exact boosts, plain addition

Monty Dabas. 8 October 2026. Python 3.12, standard library only. Exact
integer and rational arithmetic.

The owner's question: pure numbers are rotation and repetition once a
unit is fixed; is there a place for the prime numbers in the native
algebra?

Sources read before building: **EMK-1** (the block over the rationals,
R² = −1, K² = +1), **SY1** (three sectors), **PS1** (Bindu: the lift
z → z/√det z to the unit quadric), **RMG2-T1** (number sectors),
**BC1** (the boundary sector is addition), **QC1-T3** (silence),
**CL1 / R38** (clocks with Q marks; the eight-mark value 1/√2).

Carrier: rational elements z = a + bR, det z = a² + b². An *exact turn*
is a rational point of det z = 1.

## Results

**P1 — the square of the Bindu lift needs no root.**
B₂(z) = z / z̄ = z² / det z is always an exact turn.

**P2 — prime turns.** For a prime p = x² + y² (that is, p ≡ 1 mod 4)
put u_p = B₂(x + yR):

```text
u₅ = (3 + 4R)/5 ,   u₁₃ = (5 + 12R)/13 ,   u₁₇ = (15 + 8R)/17 ,   u₂₉ = (21 + 20R)/29 , …
```

Every exact turn is, in exactly one way,

```text
R^k · ∏ u_p^{n_p} ,        k ∈ {0,1,2,3},   n_p whole numbers ,
```

and the exponents add under composition. The exact turns of the native
algebra are the quarter turns times **one independent turn for each
prime p ≡ 1 mod 4**. (Checked: all 48 exact turns of denominator below
150 factor and rebuild.)

**P3 — the other primes.** A prime p ≡ 3 mod 4 is not a sum of two
squares: it gives no turn, only scale. The prime 2 = det(1 + R) gives
B₂(1 + R) = R, the quarter turn. Its half — the eighth turn — needs √2:
the eight-mark value of R38.

**P4 — no exact turn repeats.** u^n = 1 for an exact turn only if u is a
quarter turn. So the clocks that are exact over the rationals have
Q = 1, 2 or 4 marks, and no others (the cosine alone is rational also
for Q = 3 and 6). A clock that ticks by a prime turn never returns.

**P5 — split sector.** Exact boosts
h(t) = (t + 1/t)/2 + ((t − 1/t)/2)K, t rational, compose as
h(t)h(s) = h(ts): one independent boost for **every** prime.

**P6 — dual sector.** 1 + sN composes by adding s, and every element
has every root: no primes at all.

## Reading

```text
sector      exact elements compose like          primes
circular    quarter turns × whole-number vectors  one turn per prime ≡ 1 mod 4 ; 2 is the quarter turn ; primes ≡ 3 mod 4 silent
split       non-zero rationals under ×            one boost per prime, all of them
dual        rationals under +                     none
```

The owner's sentence holds in this form: once the unit is fixed
(det = 1), the exact pure numbers of the circular sector *are*
rotations, and they are built from prime rotations by repetition
(whole-number powers), uniquely. The primes are in the algebra as its
independent elementary turns.

Two consequences inside the framework:

- Every certificate in this line that used (3/5, 4/5), (5/13, 12/13),
  (−7/25, 24/25) was using u₅, u₁₃ and u₅².
- QC1's silence needs a return. By P4 no exact turn returns except the
  quarter turns; recordable cycles with other counts leave the
  rationals. The dyadic clocks of CL1 are the prime 2.

## What is not shown

- P2–P6 are the arithmetic of sums of two squares and of the rationals,
  known since long; what is the framework's own is their placement in
  the three sectors and the link to Bindu and to the clocks.
- No pure number of physics (1/137 or a mass ratio) is obtained. A
  prime turn's angle is not a rational part of a full turn, and no
  relation between those angles and measured constants was looked for
  or found.
- "Field" in the owner's question: the exact turns form a group, not a
  field; the field is the rationals with R, in which the primes ≡ 1
  mod 4 split into the two factors that make a turn.

## Claim boundary

```text
B₂(z) = z/z̄ IS AN EXACT TURN                                            PROVED
UNIQUE FACTORISATION OF EXACT TURNS INTO QUARTER AND PRIME TURNS        PROVED on all turns of denominator < 150 and on composites; classical in general
PRIMES ≡ 3 mod 4 GIVE NO TURN; 2 GIVES THE QUARTER TURN                 PROVED (to 120) / exact
NO EXACT TURN REPEATS EXCEPT QUARTER TURNS; RATIONAL COSINE ⇒ Q ∈ {1,2,3,4,6}   PROVED
SPLIT SECTOR: ALL PRIMES; DUAL SECTOR: NONE                              PROVED
A PHYSICAL PURE NUMBER FROM THE PRIME TURNS                              NOT OBTAINED
```

## Reproduce

```text
python pt1_prime_turns.py
python -m unittest test_pt1_exact
```
