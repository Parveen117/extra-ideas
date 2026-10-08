# PM1 — The prime-turn series: the series the native algebra carries in place of the draft's three prime objects

Monty Dabas. 9 October 2026. Python 3.12, standard library only. Exact integer arithmetic; directed rational
enclosures for the seam law.

Source objects, from the catalog draft:

```text
prime lattice         { (log p , 2π/p) }
prime-modular form    Σ_p e^(2πiτ/p)                    "converges for Im τ > 0", modularity conjectured
seam zeta function    Σ_p p^(−s) e^(−2πi/p)             functional equation and zeros "on the seam" conjectured
```

None of the three is kept (P6). Their idea — primes entering a series through turns — is kept, with the turns
that the native algebra has: the prime turns of PT1.

Sources read before building: **physics PT1** (exact turns; u_p = π/π̄ for each prime p ≡ 1 mod 4; 2 is the
quarter turn; primes ≡ 3 mod 4 give no turn; no exact turn repeats), **RKF F00-E** (cut-complex field),
**Publications LAM-1 T3** (theta seam at a point by directed rational enclosures), **LAM-2** (a seam phase
must be derived from finite arithmetic; the character mod 4; T4: a twist with no conductor has no such
derivation), **LAM-3** (the crossing). The native primes and Euler product are RKF F00-H / F00-I (pointed to,
not re-read for this stage).

Carrier: whole cut-complex numbers α = a + bι, norm N(α) = a² + b², units ι^j.

## Results

**P1 — among the powers, the unit-free ones are the fourth powers.** α^m is the same for the four unit
multiples of α exactly when 4 | m. For every other m each norm shell sums to zero. So

```text
A_k(n) = ¼ · Σ_(N(α) = n) α^(4k)
```

is a whole number with no turn-part (checked for n ≤ 1500, k = 1, 2).
A₁ = 1, −4, 0, 16, −14, 0, 0, −64, 81, 56, 0, 0, −238, …

**P2 — the primes enter as their turns.** A_k is multiplicative, and

```text
A_k(2) = (−4)^k                                   2 is the quarter turn :  (1 + ι)² = 2ι
A_k(p) = 0                    p ≡ 3 mod 4         no turn (PT1-P3)
A_k(p) = 2·p^(2k)·rad( u_p^(2k) )   p ≡ 1 mod 4   u_p the prime turn :  u₅ = (3 + 4ι)/5 , u₁₃ = (5 + 12ι)/13 , …
A_k(p^(e+1)) = A_k(p)·A_k(p^e) − χ(p)·p^(4k)·A_k(p^(e−1))        χ = the character mod 4 .
```

**P3 — size.** |A_k(p)| < 2p^(2k): the rad-part of a turn is below 1, and it is never ±1 because no prime
turn repeats (PT1-P4).

**P4 — product form.** The numbers built from the prime rule P2 alone equal the shell sums for every
n ≤ 1500 (two routes). Formally

```text
Σ A_k(n)·n^(−s) = ∏_p 1 / ( 1 − A_k(p)·p^(−s) + χ(p)·p^(4k − 2s) ) .
```

**P5 — seam law at a point.** Θ_k(t) = Σ_α α^(4k)·Exp(−π·N(α)·t), which for k ≥ 1 is 4·Σ A_k(n)·Exp(−πnt)
(for k = 0 the term α = 0 adds 1). The law is

```text
Θ_k(1/t) = t^(4k+1) · Θ_k(t) ,        phase +1 .
```

At t = 2 and t = 3/2, for k = 0, 1, 2, the two sides lie in two-sided rational enclosures that overlap and
are narrower than 10⁻⁴⁰ (π by two arctangent series, Exp by the factorial series, all rounding directed).
That certifies agreement to forty digits at these points; equality is the classical anchor, pinned.
Θ₁(½) = 0.2372471752… ; k = 0 is the square of LAM-1's theta. Separated by the same enclosures: the weights
t^(4k) and t^(4k+2), and the phase −1.

Where the phase comes from. The anchor (classical, pinned) is: a harmonic reading of degree d crosses with
weight t^(d+1) and a phase that is the d-th power of a quarter turn. The harmonic readings of degree d are the
rad- and turn-parts of α^d; by P1 they survive the sum over a shell only for d = 4k, and then the phase is 1.
For the other degrees the phase is not 1 and the series is identically zero: the two statements agree.
Being unit-free is not enough by itself: a⁴ + b⁴ = ¾N² + ¼·rad(α⁴) is unit-free, is not a power, and does
not cross (gap 1.11).

**P6 — why the draft's series are not used.**

```text
Σ_p e^(2πiτ/p)            |term| = Exp(−2π·Im τ/p) > ½ for p > 14·Im τ : the terms do not tend to zero ;
                          the series converges for no τ. (The modulus the draft quotes tends to 1, not to 0.)
twist e^(−2πi/p)          a turn of finite order p ; no exact turn other than a quarter turn repeats (PT1-P4) ;
                          as a turn 1/n it is not multiplicative (½ + ⅓ ≠ ⅙) ; it differs at every prime, so it
                          depends on no residue class : no conductor (compare LAM-2 T4).
Σ_p p^(−s) e^(−2πi/p)     differs from the plain prime sum Σ_p p^(−s) by a series that converges absolutely for
                          Re s > 0 (|Exp(−ιθ) − 1| ≤ θ) : the twist adds nothing singular there. The plain prime
                          sum has logarithmic branch points inside 0 < Re s ≤ 1 (at s = 1, at 1/k, and at ρ/k for
                          the zeros ρ) and cannot be continued past Re s = 0 (classical ; recalled, not certified
                          here). No s ↔ 1 − s law is known or derived for either sum.
(log p , 2π/p)            the native pair is (p , u_p), PT1 ; it is what P2 reads.
```

The draft's expansion Σ_p p^(−s)e^(−2πi/p) = P(s) − 2πi·P(s+1) − 2π²·P(s+2) + … is correct for Re s > 1.

## Use

- It is the next case of LAM-2's rule (a seam phase must be derivable as finite seam arithmetic): after the
  characters mod q, the characters of the cut-complex whole numbers themselves. The data per prime are an exact
  turn, the product form is exact, and the phase of the crossing is fixed by the four units.
- One object carries both the product over primes (P4) and a seam law (P5). The draft's seam zeta has neither.

## What is not shown

- P2–P5 are known mathematics: the theta series of the square lattice with a harmonic coefficient and the
  L-series of a character of the Gaussian whole numbers (Hecke; recalled, not re-read at source). The line's
  part: the reading through PT1's prime turns, the four units fixing the degree and hence the phase, the
  certificates, and P6.
- P5 is agreement at the stated points, not a proof for every t and not a proof of equality. Its anchor is
  classical and pinned, as in LAM-1/2. It is a property of the square lattice with a harmonic unit-free
  reading; the primes enter through P2 and P4.
- Multiplicativity is checked to n = 1500; in general it rests on unique factorisation of cut-complex whole
  numbers (classical; PT1 states the same boundary).
- No continuation of Σ A_k(n) n^(−s) beyond its region of convergence is made here, no completed function,
  no statement about zeros. N1, N2, N3 and RH are untouched and OPEN.

## Claim boundary

```text
AMONG THE POWERS, UNIT-FREE = FOURTH POWERS ; A_k WHOLE, NO TURN-PART                                PROVED
A_k MULTIPLICATIVE ; PRIME VALUES = PRIME TURNS ; PRIME-POWER RULE WITH χ MOD 4 ; PRODUCT FORM        PROVED to n = 1500 ; classical in general
|A_k(p)| < 2p^(2k)                                                                                   PROVED
Θ_k(1/t) = t^(4k+1)·Θ_k(t) AT t = 2, 3/2 ; k = 0, 1, 2 ; PHASE +1                                    AGREEMENT TO 40 DIGITS CERTIFIED (enclosures) ; equality classical, pinned
THE DRAFT'S PRIME-MODULAR FORM (DIVERGES) AND SEAM ZETA (NO SEAM LAW)                                NOT KEPT
CONTINUATION, COMPLETED FUNCTION, ZEROS ; N1 / N2 / N3 ; RH                                          NOT CLAIMED ; OPEN
```

## Reproduce

```text
python pm1_prime_turn_series.py
python -m unittest test_pm1
```
