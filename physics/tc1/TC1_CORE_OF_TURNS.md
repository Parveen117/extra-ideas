# TC1 — The core of d turns: three numbers, a difference form whose diagonal is the valley, a flat measure at four turns, the constant of the logarithm, and the scale the core makes for itself

Monty Dabas. 9 October 2026. Python 3.12, standard library only. Exact rationals throughout.

TV1 found that the quadratic layer of the constant modes closes for three turns and leaves a logarithm at the
centre points for four. The question put to this stage: at which scale does the quartic core cut that logarithm
off, and does the core fix the size of a mass or only its form.

Sources read before building: **TV1**, **HL1** (H7), **TW1**; physics **SC1** (scale is the unit of count; no
law in the block gives a dimensional constant), **DM1** (the reading of the missing cut), **DS1** (the count
past the diagonal); RKF **theorum/24** (F = R − D), **theorum/28** (outward form); **YM99** and the record of
YM100 (non-quadratic core of the constant modes, scale 1/L).

## The record

Near a centre point the vector parts of d unit blocks are a free matrix C of three rows (the cut components)
and d columns (the directions). The commutator weight of all pairs is X = Σ_(i<j) |c_i × c_j|² (TV1-V1); the
core is this weight with nothing along a common axis removed. G = C·Cᵗ is the 3 × 3 matrix of the rows.

## Results

**C1 — three numbers.**

```text
X = ½·[ (tr G)² − tr G² ] = x₁x₂ + x₂x₃ + x₃x₁ ,        x₁, x₂, x₃ the squared singular values of C .
```

The weight is the second symmetric function of three numbers, for any number of directions. It is unchanged by
turning the three cuts and by turning the d directions (exact rotations). A fourth direction adds no fourth
number: every 4 × 4 minor of CᵗC is zero. Three cuts read three directions at once; more directions only turn
into one another. (100 exact matrices, d = 2 … 6.)

**C2 — it is the one law, and the valley is its diagonal.**

```text
X = R − D ,      R = ½(x₁ + x₂ + x₃)² ,      D = ½(x₁² + x₂² + x₃²) ,      0 ≤ X ,      ⅓ ≤ D/R ≤ 1 .
```

X is zero exactly where D = R: one number alone, the valley. Lost over seen is 1 there and ⅓ where the three
numbers are equal; it is never above 1, so the difference form never changes sign on the record — the count
past the diagonal (DS1) is zero. As a quadratic form in (x₁, x₂, x₃) it has values 2, −1, −1: signature (1, 2),
and the valleys are its null lines.

**C3 — the layer about a valley point.** With C = **n**·pᵗ + T (T across the axis **n**, p the amplitudes of
the d turns along it):

```text
X = |p|²·|T|² − |T·p|² + X(T) .
```

In the space of turns the layer is |p|² times the projector across p: every transverse direction has the same
rate, and the one direction along p is free (turning the axis). The ground rate of one transverse direction is
|p|/√2, so the layer rises along the valley as √2·(d − 1)·|p| — linearly. In TV1's variables its determinant
is the tree determinant (Π p_i²)(Σ p_i²)^(d−2).

**C4 — the number of directions is one exponent, and four is where it vanishes.** In the three numbers the
free measure of C is

```text
Π_(a<b) |x_a − x_b| · Π_a x_a^((d−4)/2) · d³x .
```

For four turns the squared singular values are flat coordinates. Checked on the means of e₁, e₂, e₃, e₁²,
e₁e₂ under the free weight, from the entries of C (pairing) and from this measure: 6, 9, 3, 42, 72 at d = 4,
and the same agreement at d = 6 with the factor x₁x₂x₃; the flat measure does not give the means of six.

**C5 — the constant of the logarithm.** Four turns, flat core, the largest number cut at Λ:
K(Λ) = ∫_(Λ > x₁ > x₂ > x₃ > 0) Π(x_a − x_b)·Exp(−2X)·d³x. Along the valley (x₁ large, x₂, x₃ of order 1/x₁)
the inner integral is at most 1/(16·x₁), and at least 1/(16x₁) − 9/(64x₁³) + … (the weight's own correction,
compared term by term). So

```text
(Log 2)/16 − δ(Λ)  ≤  K(2Λ) − K(Λ)  ≤  (Log 2)/16 ,        0 < δ(Λ) < 1/(8Λ²) .

Λ          2          4          8          16         32
δ(Λ)       0.0181     0.00317    0.00082    0.00021    0.00005
```

Every doubling of the largest number adds (Log 2)/16: the record has density dx₁/(16·x₁) along the null line.
The measure has degree 6 in the numbers and X degree 2, so with a coupling κ in the weight the record is
κ⁻³·K(Λ·√κ): the core's own scale is x ~ κ^(−1/2), that is |C| ~ κ^(−1/4), and the coupling enters the
logarithm with half the weight of the size — (Log κ)/32 against (Log Λ)/16.

**C6 — the scale the core makes for itself.** X has degree four in C. For the Hamiltonian of the constant modes
in the core, −(g²/2)·Δ + X/g², the change C → g^(2/3)·C gives

```text
−(g²/2)·Δ + X/g²  =  g^(2/3) · ( −½·Δ + X ) .
```

Every rate of the core is g^(2/3) times a pure number; no scale is put in. The layer's ground rate |p|/√2 does
not depend on g at all: the valley is closed by a rise that carries no coupling. A Gaussian reading gives the
ground rate at most c_d·g^(2/3) with c_d³ = (729/256)·d³(d − 1): c₃ = 5.36, c₄ = 8.18 (upper bounds).

## What this answers

```text
at which scale is the logarithm cut off          at the core's own: |C| ~ κ^(−1/4) in the record, g^(2/3) in the rates ;
                                                 the scale follows from the degree of the commutator, nothing is put in
does the cut fix the size of a mass              no. It fixes the form: the power 2/3 and pure numbers. The size is the
                                                 cell's: a rate of the constant modes is (pure number)·g^(2/3) in units of
                                                 the torus (YM100's 1/L). SC1 said the same of every law in the block.
which four, which three                          the block has three cuts, so any number of directions is read through
                                                 three numbers; the number of directions is the exponent (d − 4)/2 of
                                                 the measure, zero at four
```

## What is not shown

- No mass gap. The pure numbers of −½Δ + X (its lowest rates, their difference) are not computed; only an upper
  bound on the ground rate is given. A lower bound on the difference needs an outward certificate (theorum/28)
  that is not built.
- A rate proportional to g^(2/3) over the size of the torus goes to zero with the volume at fixed coupling. That
  the coupling of the constant modes changes with the size, and how the logarithm of C5 enters there, is not
  shown. This is the open wall of MG1's verdict, now with its form and constants.
- C5 is the flat core with a sharp cut at Λ. For the compact record of four unit blocks the cut is the block
  itself; the constant there is not certified.
- C4 is checked on five means at d = 4 and d = 6, not proved as a change of variables.
- In older terms (general knowledge, not re-read): C1 is the sum of squared 2 × 2 minors; C4 is the eigenvalue
  measure of a matrix of squares; the g^(2/3) law and a linearly rising valley are known for gauge fields in a
  small box. The line's part: the weight as the one law with the valley as its diagonal and lost/seen between ⅓
  and 1; that a fourth direction adds no number and the measure is flat exactly at four; the constant (Log 2)/16
  with certified bounds; the layer as one projector in the space of turns with a rise free of the coupling; and
  the answer "form, not size".

## Claim boundary

```text
X = e₂ OF THREE NUMBERS ; INVARIANT UNDER BOTH TURNINGS ; NO FOURTH NUMBER            PROVED (exact, d ≤ 6)
X = R − D ≥ 0 ; ZERO ON THE VALLEY ; ⅓ ≤ D/R ≤ 1 ; SIGNATURE (1, 2)                    PROVED (exact)
LAYER = |p|²·(PROJECTOR ACROSS p) ; RISE √2(d − 1)|p| ; TREE DETERMINANT              PROVED (exact)
MEASURE Π|x_a − x_b|·Πx^((d−4)/2) ; FLAT AT FOUR                                      CHECKED on five means, d = 4, 6
(Log 2)/16 − δ ≤ K(2Λ) − K(Λ) ≤ (Log 2)/16 ; κ⁻³K(Λ√κ)                                PROVED (flat core, exact bounds, Λ = 2 … 32)
RATES = g^(2/3) × PURE NUMBERS ; LAYER'S RISE FREE OF g ; c₃ < 5.36, c₄ < 8.18         PROVED (the bounds are upper bounds)
THE PURE NUMBERS ; A GAP UNIFORM IN THE VOLUME ; FOUR DIMENSIONS                      NOT COMPUTED ; NOT SHOWN ; OPEN
```

## Reproduce

```text
python tc1_core_of_turns.py
python -m unittest test_tc1
```

## Later note (CR1, 9 October)

[CR1](../cr1/CR1_CORE_RATES.md) computes the pure numbers left open in C6, in the sector unchanged by both
turnings: lowest rate in [4.14, 5.1868] and second below 8.0480 for three turns; [6.32, 8.0030] and 11.5660 for
four (exact counts on a ladder; an elementary bound from below). The difference, 2.86 and 3.56, is not a
certified bound. Nothing above is changed.
