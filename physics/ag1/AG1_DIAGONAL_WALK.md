# AG1 — The Riemann line's seam on the Yang–Mills record; the diagonal walk is the arithmetic–geometric mean

Monty Dabas. 9 October 2026. Python 3.12, standard library only. Exact integer series; directed rational
enclosures (tools/pm1).

The owner's instruction: carry the Riemann line's results to Yang–Mills — the same along the seam on the real
and on the imaginary axis, and after that the diagonal walk. If it closes, good; if not, the direction.

Sources read before building: tools **DS1** (F = Z*𝔧Z; the diagonal vector 1 + ι is the prime 2; the mirror
form; no crossing), **PM1** (the lattice of whole cut-complex numbers and its seam law), **RW1**; Publications
**LAM-1 T3, LAM-2, LAM-3** (theta seam by enclosures; the character mod 4; cancellation crossed by exact
arithmetic), **YM75** (the Riemann line's Schur and Birman–Schwinger tools already carried to the actual
Yang–Mills gap — not repeated here); physics **DU1** (doubling the cell; the record that is its own dual),
**CT1**, **CZ1** (the weight of a twist; weak-end exponent π²/2F), **MG1** (the weak-coupling wall); RKF
**theorum/24**, **theorum/28**.

## The record

A pair of turns with the heat weight. Its contents are the whole cut-complex numbers α = a + bι — the two axes
are the two turns. A content weighs q^N(α), q = Exp(−πt), t the heat time of the cell (small t: weak coupling).

```text
S(t) = Σ q^(a²+b²)                       source
F(t) = Σ (−1)^(a+b) q^(a²+b²)            flip: contents on the diagonal sublattice (1 + ι)·ℤ[ι] minus the rest
L(t) = Σ q^((a+½)²+(b+½)²)               the half-shifted sum
```

F is theorum/24's R − D for that cut; it is the centre mark of the record (the sign of a + b).

## Results

**A1 — the seam, the same on both axes.** Under the mirror t → 1/t (contents ↔ windings; the Riemann line's
theta seam, LAM-1 T3 and PM1-P5 at k = 0):

```text
S(1/t) = t·S(t) ,        F(1/t) = t·L(t) ,        L(1/t) = t·F(t) .
```

At t = 1 the lattice of contents and the lattice of windings are the same square lattice: real axis and
imaginary axis alike. Enclosures at t = 2, 3/2, 5, narrower than 10⁻¹⁰⁰; a wrong partner is separated.

**A2 — a right triangle.**

```text
S² = F² + L²        (exact, as series of whole numbers) .
```

The mirror exchanges the two legs. On the seam they are equal: F = L and (F/S)² = ½ — the eight-mark value
1/√2, an angle of 45°. Off the seam they are not (t = 1.1: L < F).

**A3 — the diagonal step.** Multiplying the lattice by the diagonal vector 1 + ι doubles every norm and lands on
the even coset. So

```text
R(q) = (S + F)/2 = S(q²)         the seen part of the record is the record of the doubled cell ,
count(norm 2n) = count(norm n) .
```

The mirror J (exchange of a and b) cancels the odd coset in the mixed flip Σ(−1)^b q^N, which is therefore the
flip of the doubled cell; its square is S·F.

**A4 — the walk is the arithmetic–geometric mean.** For one doubling of the cell (t → 2t):

```text
S′ = (S + F)/2 ,        F′ = √(S·F) .
```

Hence, at every cell size and every coupling,

```text
F < F′ < S′ < S                          the flip only grows, the source only falls, they never cross ;
S′² − F′² = ((S − F)/2)²                 the doubled cell's source is the hypotenuse over its flip and the lost part ;
M(S, F) = 1                              the arithmetic–geometric mean of source and flip is the same at every
                                         cell size: it is 1, the count of the origin .
```

From t = 1/8 (F/S = 1.4·10⁻⁵) nine steps of the mean give 1 to fifty digits. The seen and the lost change at
every step; this one number of the pair does not.

**A5 — the weak end is the diagonal, and the mirror crosses the wall.** As t → 0, F/S → 0: seen = lost. The
flip there is

```text
F(t) = (4/t)·Exp(−π/(2t))·(1 + Exp(−2π/t) + …)² ,
```

smaller than every power of t (at t = 1/60 it is below t²⁰). Described by contents it is what is left of a sum
of terms of size 1: 11 digits cancel at t = 1/20 and 38 at t = 1/60. Described by windings it is one positive
term. Both descriptions are enclosed and agree; the first winding term is the flip to a thousandth.

```text
t = 1/20        F = 1.816881·10⁻¹²
t = 1/60        F = 2.811387·10⁻³⁹
```

**A6 — the turn block.** Contents m = 2j + 1 with weight m² and heat weight (u = t/4π). A closed surface has
flip K₋(u)/K₊(u), K₊ = Σ m²·Exp(−πu·m²), K₋ = Σ(−1)^(m−1) m²·Exp(−πu·m²). In the mirror description

```text
K₊ = 1/(4π u^(3/2)) · Σ_n (1 − 2πn²/u)·Exp(−πn²/u) ,
K₋ = 1/(4π u^(3/2)) · Σ_n (2π(n+½)²/u − 1)·Exp(−π(n+½)²/u) .
```

The two descriptions agree (u = 1, 1/4, 1/16, 1/36); the flip is positive and below 1. At the weak end it is
the first winding term, 2(π/2u − 1)·Exp(−π/4u): 1.68·10⁻⁴ at u = 1/16, 5.84·10⁻¹¹ at u = 1/36. Its exponent is
π²/T in the heat time T of the whole surface. Matching the heat weight to the Wilson weight near the identity
(Exp(−x²/t) against Exp(−κx²/2)) gives t = 2/κ for a face and T = 2F/κ for F faces: the exponent is π²κ/(2F),
the weak-end exponent CZ1 found for the Wilson weight.

**A8 — the turn block walks at least as fast.** Let ρ = K₋/K₊ be the flip of the closed surface at heat time T
and ρ′ the same at 2T. The pair record has ρ′ = 2√ρ/(1 + ρ) exactly (A4). The turn block has

```text
ρ′·(1 + ρ)  >  2√ρ
```

at every coupling tried: strictly, by enclosures, at nine couplings from u = 1/100 to u = 2. The ratio of ρ′ to
the mean step is 4.39, 3.71, 2.59, 1.67, 1.15, 1.003 at u = 1/100 … 1/4 and tends to 1 at the strong end, where
the difference begins with 64·q¹² (q = Exp(−πu); exact series). The non-abelian block leaves the diagonal at
least as fast as the pair. This is checked, not proved for all couplings.

**A7 — the reflection form is the mirror form.** For a transfer record c(m) = Σ w·λ^m, the form
Σ f_s f_t† c(s+t) has inertia (real values + pairs, pairs): DS1-D4 with the mirror λ → λ†. It is non-negative
exactly when every transfer value is on the rad axis. The gap is the same count with the line moved:
n₋(r·c(s+t) − c(s+t+1)) = the number of values above r.

## The two lines side by side

```text
Riemann line                                   this record
mirror s → 1 − s† ; line rad = ½               mirror t → 1/t ; seam t = 1                       A1
(the two are Mellin partners: the seam law S(1/t) = t·S(t) is the law s ↔ 1 − s of the source's series)
Euler factor at the prime 2                    one diagonal step: q → q², multiplication by 1 + ι   A3
n^(−½): the diagonal share                     the diagonal sublattice is half the lattice          A3
the alternating series                         the flip F                                            A3
positivity of a mirror form ⇔ all on the line   reflection form ≥ 0 ⇔ all transfer values real      A7, DS1-D4
count k_Σ = N₊(B − 1)                          gap count n₋(r·c − c′)                                A7, DS1-D3
never across (to be shown for the five-matrix) never across: F < S at every cell and coupling       A4
```

The series of the source is the count of whole cut-complex numbers of each norm, 4·Σ_(d|n) χ(d) with χ the
character mod 4 (checked to n = 400): the product of the Riemann series and LAM-2's series. The two lines are
here two readings of one record.

One difference of kind (A7, DS1-D4). On the Riemann line the open statement is *across* the seam (every zero
on it); the statement *along* it (the first zero lies at a positive height) is known. On the Yang–Mills line the across statement — transfer values on the rad
axis, reflection positivity — holds by construction, and the open statement is *along*: the first value stays
a finite distance from 1 as the cell shrinks. The Riemann line's count tools therefore carry to the
Yang–Mills line where positivity must survive a limit (YM55, YM75), and the gap needs the walk.

## What it gives and what it does not

```text
for the pair record            the flip is positive at every coupling; the walk is exact (A4) ; the weak end is the
                               diagonal, reached only in the limit ; the flip there is 4/t·Exp(−π/2t), with no
                               expansion in t ; the wall of the content description is crossed by the mirror      PROVED
for a closed surface of the    flip positive at every coupling, weak end = first winding term (A6) ;              PROVED (heat weight)
turn block (two dimensions)    its walk is at least the mean step of the pair (A8)                                CHECKED at nine couplings
for the fabric in four         NOT SHOWN. There is no exact block step. See the direction below.
dimensions
```

## Direction (a proposal, not a result)

A4 is two statements about one block step: the source goes to the arithmetic mean, the flip to the geometric
mean. Never-across is then AM ≥ GM, and the invariant is their common mean. For the four-dimensional fabric the
same conclusion needs only two inequalities for one block step — source at most the arithmetic mean, flip at
least the geometric mean, each with a controlled remainder — with theorum/28's outward form for the remainder.
The exact pieces at hand: at fixed cosets the dual couplings add (CZ1-Z3); every face has the mirror
description (A6), which is the convergent one at the weak end; YM75 keeps the whole source matrix; and the
one non-abelian record with exact data obeys the flip inequality with room (A8). Decimation
inequalities of this kind have been proposed for this problem in the literature and are not accepted as
complete (general knowledge; not re-read). What this stage adds to that direction is the exact model in which
both hold with equality, the invariant, and the reading of a small mass count: a flip of size Exp(−c/t) at the
weak end, which no expansion around the diagonal can see and the mirror sees as one term.

## What is not shown

- No statement about the mass gap of the four-dimensional theory, and none about the Riemann hypothesis.
- A1 and A6 are agreement of enclosures at the stated couplings; equality for all t is the classical seam law,
  pinned as in LAM-1. A2–A4 are exact identities of series (checked to order 400) with the written lattice
  proofs above; they are Jacobi's and Landen's identities and Gauss's mean (known since long). The line's
  part: the reading as source, flip and lost part of one cut, the diagonal vector as the block step, the
  never-across and the invariant in the law's terms, and the match with CZ1 and the Riemann line.
- "Mellin partners" in the table is the classical relation between a seam law in t and a law s ↔ 1 − s; it is
  stated, not certified here.
- The pair record is abelian. That the four-dimensional non-abelian fabric obeys the two inequalities is not
  known.
- The heat weight is used throughout. CZ1's Wilson weight has the same weak-end exponent; no more is claimed.

## Claim boundary

```text
S(1/t) = t·S ; F(1/t) = t·L ; L(1/t) = t·F                                   AGREEMENT CERTIFIED at t = 2, 3/2, 5 ; equality classical, pinned
S² = F² + L² ; ON THE SEAM F = L, (F/S)² = ½                                  PROVED (series) ; seam value to 100 digits
R(q) = S(q²) ; WALK S′ = (S+F)/2, F′² = S·F ; F < F′ < S′ < S ; M(S, F) = 1    PROVED (series ; points ; nine steps of the mean)
WEAK END: F = (4/t)·Exp(−π/2t)(1 + …)² ; TWO DESCRIPTIONS AGREE               PROVED at t = 1/20, 1/60
TURN BLOCK: MIRROR DESCRIPTION ; FLIP POSITIVE ; EXPONENT π²κ/(2F)            PROVED at u = 1, 1/4, 1/16, 1/36 (heat weight)
TURN BLOCK: ρ′(1 + ρ) > 2√ρ                                                   CHECKED at nine couplings (enclosures) ; strong-end series ; not proved in general
REFLECTION FORM = MIRROR FORM ; GAP AS A COUNT                                PROVED (finite)
THE WALK IN FOUR DIMENSIONS ; A MASS GAP ; RH                                 NOT SHOWN ; OPEN
```

## Reproduce

```text
python ag1_diagonal_walk.py
python -m unittest test_ag1
```

## Later note (TW1, 9 October)

A8 is proved at every coupling in [TW1](../tw1/TW1_TURN_BLOCK_WALK.md). The step of the turn block is exact
there — a boost by the pair's ratio, ρ′r′ = (ρ − r)/(1 − rρ) — and the record is the pair's lost parts summed
along the walk. The flip of the turn block is positive at every coupling by the walk from the strong side, so
A6's positivity no longer rests on the couplings tried. Nothing above is changed.

## Later note (SL1, 9 October)

[SL1](../sl1/SL1_SEAM_LAW_FROM_THE_WALK.md) proves the seam law of the pair — S(1/u) = u·S, F(1/u) = u·L, L(1/u) = u·F — from the diagonal
step, the two means and one count at the weak end; the only input from outside the walk is the area under
Exp(−πx²). So A1 is proved, A2 follows from the step's substitutions, and the mirror description of A6 is one change in u of that law. Nothing above is changed.
