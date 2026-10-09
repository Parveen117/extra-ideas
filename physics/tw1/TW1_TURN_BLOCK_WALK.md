# TW1 — The turn block's walk: one diagonal step is a boost by the pair's flip, the record is the sum of the lost parts, and the flip inequality holds at every coupling

Monty Dabas. 9 October 2026. Python 3.12, standard library only. Exact integer series to order 600; exact
rational sums with bounded rests; directed rational enclosures (tools/pm1, AG1).

The owner's instruction: go on from where the Yang–Mills development stood, and take the step AG1 left
open — its A8 was checked at nine couplings and not proved.

Sources read before building: **AG1** (the pair's seam and walk; the turn block A6; A8), **CZ1** (the weight of
a twist of a closed surface; d^χ), **CT1** (centre coupling, odd over even), **DU1** (dual coupling), **MG1**
(verdict: the open wall is the gap at weak coupling, YM98 §5 and YM99), **CG1** (status of theorum/41), tools
**DS1**, **RW1**, **PM1**; RKF **theorum/24** (F = R − D), **theorum/28** (outward certificate).

## The record

A closed surface of the turn block with the heat weight (AG1-A6), q = Exp(−πu), contents m = 1, 2, 3, …:

```text
K₊ = Σ m²·q^(m²)                    source
K₋ = Σ (−1)^(m−1)·m²·q^(m²)         flip                     ρ = K₋/K₊
```

and one turn of AG1's pair at the same q:

```text
a = Σ q^(n²) ,   b = Σ (−1)ⁿ q^(n²)   (n whole) ;      S = a² , F = b² ;      r = b/a ;      L² = S² − F² the lost part .
```

Term by term K₊ = −(1/2π)·da/du and K₋ = (1/2π)·db/du: the closed surface of the turn block is one turn of
the pair read through its first change in the heat time. A prime marks the next diagonal step, q → q²; L_n is
the lost part after n steps.

## Results

**W1 — one diagonal step is exact, and it is a boost.**

```text
4·a′·K₊′ = a·K₊ − b·K₋ ,          4·b′·K₋′ = a·K₋ − b·K₊ .
```

Proof on the lattice: a·K₊ − b·K₋ is the sum of m²·q^(m²+n²) over the pairs with m + n even. Those are the
pairs (x + y, x − y) — AG1-A3's diagonal sublattice — and m² = x² + y² + 2xy, where the cross term cancels
under x → −x. What is left is Σ(x² + y²)·q^(2(x²+y²)) = 4·a′·K₊′. The second line is the same with the mark
(−1)^m = (−1)^(x+y). Both are checked as series of whole numbers to order 600, and against the lattice sum.

(a′K₊′, b′K₋′) is ¼·a times the matrix [[1, −r], [−r, 1]] applied to (K₊, K₋): up to scale, a boost whose
velocity is the pair's ratio r. In ratios

```text
ρ′·r′ = (ρ − r)/(1 − r·ρ) .
```

With ρ = tanh α, r = tanh γ: tanh α′·tanh γ′ = tanh(α − γ). Rapidities subtract. (The pair itself has
tanh²γ′ = tanh 2γ, AG1-A4.)

**W2 — the record is the sum of the lost parts along the walk.**

```text
16·K₋ = b·( L₀² + 2·L₁² + 4·L₂² + 8·L₃² + … ) ,
16·K₊ = a·( L₀² − 2·L₁² − 4·L₂² − 8·L₃² − … ) .
```

So the flip of the turn block over the flip of one turn of the pair is

```text
ρ / r = (1 + Λ)/(1 − Λ) ,        Λ = (ρ − r)/(ρ + r) = Σ_(n≥1) 2ⁿ·L_n² / L₀² .
```

What the pair loses at each later step of the walk, counted 2ⁿ times and measured by what it lost at the
start. Every L_n² is a series of one sign (it begins 16·q^(2ⁿ)), so every term is positive.

This is the one law (theorum/24) with nothing added:

```text
R = L₀²   seen ,      D = 2·L₁² + 4·L₂² + …   memory ,      16·K₋/b = R + D ,      16·K₊/a = R − D ,      lost/seen = Λ .
```

Λ is below 1 at every coupling (W3), so the difference form never changes sign: the record is never across the
diagonal, and a finite number of terms with a bounded rest certifies it at any coupling — theorum/28's outward
form, here with a scalar. The diagonal seen = lost is the weak end, approached and not reached:

```text
u            1         1/4         1/16          1/64
Λ            0.0861    0.6817      0.92042253    0.98010563
1 − Λ        0.9139    0.31826     0.07957747    0.01989437
4u/π         1.2732    0.31831     0.07957747    0.01989437
```

At the weak end the distance from the diagonal is 4u/π and nothing more that a power of u can see (equal to
sixteen digits at u = 1/16).

Proof from W1 alone. (i) P = K₋/b + K₊/a obeys P′/L₁² = P/L₀² — an invariant of the walk; at the strong end
P = 2q + …, L² = 16q + …, so P = L²/8 at every coupling. Since P is the change of Log r up to 2π, this says
d Log r/du = (π/4)·L²: the pair's ratio rises with the heat time at the rate of its lost part. (ii) W1 in the variable Λ is affine:
Λ = (S − F)/(2(S + F))·(1 + Λ′), and (S − F)/(2(S + F)) = 2·L₁²/L₀² by AG1-A4. Each factor is below ½ and
0 < Λ < 1 (W3), so the sum converges with a rest below its last term. All four lines are checked as series to
order 600; at u = 1/4 twelve lost parts give ρ/r to 60 digits (Λ = 0.6817; the first four terms carry 72.3,
25.5, 2.2 and 0.008 per cent).

**W3 — positivity is carried by the walk from the strong side.** For q ≤ ½ the terms of b and of K₋ fall, so
b > 1 − 2q ≥ 0 and K₋ > q − 4q⁴ > 0. One step back keeps both: b·a = b′² (AG1-A4) and
a·K₋ = 4·b′·K₋′ + b·K₊ (W1). Every q below 1 reaches q ≤ ½ after finitely many steps. Hence, at every
coupling and without the mirror description,

```text
0 < r < ρ < 1 ,        and        3γ < α < γ + γ′ ,
```

that is (3r + r³)/(1 + 3r²) < ρ < (r + r′)/(1 + r·r′): the first from ρ′ > r′, the second from ρ′ < 1. The
lower bound is above 2r/(1 + r²) = F′/S′ (the difference is r(1 − r²)² over a positive number): the flip of
the turn block at q is larger than the pair's flip at q².

```text
u        r          ρ             lower bound    upper bound
1/4      0.0864     0.4567        0.2542         0.4833
1/2      0.4142     0.9306        0.8673         0.9309
1        0.8409     0.99935461    0.99870967     0.99935462
```

**W4 — the weak end by walking back.** With β = ρ/r, W1 read backward is

```text
β + 1 = 2·(β′ + 1)/(1 + δ) ,        δ = 2r²·β′/(1 + r²) ,        r(√q) ≤ r(q)² .
```

Each step back doubles β + 1 up to a factor 1 + δ, and δ falls doubly exponentially because r at least squares.
So the pair's flip and the turn block's flip are both smaller than every power at the weak end, while their
ratio grows only as a power: (β + 1)/2ⁿ has a positive limit. Walked back from u = 1 with enclosures:

```text
the flip at u = 1/16 :  1.6832·10⁻⁴ , equal to AG1-A6's mirror description to 80 digits ;
the limit of u·(β + 1) :  π/2 to 80 digits — the mirror's first winding term .
```

The description by windings is not used to get there. In W2's terms 1 − Λ = 2/(β + 1): at the weak end the
later lost parts add up to the first one, and the source of the turn block is what is left over.

**W5 — the flip inequality at every coupling.**

```text
ρ′·(1 + ρ)  >  2·√ρ        for every u > 0 .
```

The pair has equality (its step is the geometric over the arithmetic mean, AG1-A4). Proof in two parts that
overlap.

*Strong side, q ≤ 7/10.* With K₊ + K₋ = 2·O (odd m), K₊ − K₋ = 8·K₊(q⁴) the inequality is

```text
k(q²)²·k(q⁴)²  >  o(q)²·o(q²)·k(q⁸) ,        k(x) = K₊(x)/x = 1 + 4x³ + 9x⁸ + … ,   o(x) = O(x)/x = 1 + 9x⁸ + 25x²⁴ + …
```

(an identity of series, order 600). The difference of the two sides begins 8q⁶ − 18q⁸ and is not of one sign
term by term (73 negative coefficients to order 600), so there is no proof by the sign of its terms. But each
side is a series of one sign. Left − 1 ≥ 8q⁶. Right − 1 is q⁸ times a rising series, and at q = 3/5 it is below 8q⁶:
that settles every q ≤ 3/5. From 3/5 to 7/10, in forty pieces of width 1/400, the left side at the lower end
exceeds the right side at the upper end (smallest margin 2.8 per cent). Exact rationals; sixteen contents and a
bounded rest.

*Weak side, q ≥ 7/10.* In (r, β) the inequality reads (1 + r²)(β − 1)²(1 + rβ)² > 8β(1 − r²β)², which holds
whenever (β − 1)² > 8β, so for β ≥ 99/10. Every q ≥ 7/10 is reached by at least two steps back from the base
[(7/10)⁴, (7/10)²]. On each of 100 boxes of the base, r and β are enclosed by exact sums; by W4,
β + 1 ≥ 4(β₀ + 1)/Π(1 + δ) with the product bounded once and for all. The certified bound is β ≥ 12.67 for
every q ≥ 7/10 (the value at 7/10 is 12.84).

```text
u                          0.01     0.05     0.1      0.2       0.3        1
ρ′(1 + ρ)/(2√ρ)            4.39     1.89     1.28     1.016     1.0005     1 + 8q⁶ − …       (floating, for the eye)
```

**W6 — what follows.** The mean step is above ρ, so

```text
ρ(q²) > ρ(q)                               the flip rises at every step, at every coupling ;
ε(q²) < ε(q)² ,   ε = (1 − ρ)/(1 + ρ)      centre-odd over centre-even (CT1's tanh k_c of the closed surface)
                                           falls below its square at every step .
```

With the heat weight a closed surface of sphere type with N faces at face time t is this record at u = N·t
(CZ1-Z2 with χ = 2; AG1-A6). So for a surface with twice the faces the odd over even is below the square:
ε(2ⁿ·u) < ε(u)^(2ⁿ). An exponential law
in the area whose rate −Log ε(u)/u can be read at any size and at any coupling, and only improves under the
walk.

```text
u          1/8       1/4       1/2        1
ε          0.9172    0.3730    0.0359     3.23·10⁻⁴
ε at 2u    0.3730    0.0359    3.23·10⁻⁴  2.60·10⁻⁸
```

Controls. For the torus record (weights 1 in place of m²; χ = 0 in CZ1's d^χ) the inequality is reversed at
q = 1/5, 1/3, 1/2 (exact enclosures): the turn block's m² is needed. A claim five per cent stronger is not
certified by the grid; a bound of 12.87 is not certified on the weak side; the step with the wrong sign fails as
a series.

One number for the eye (floating, not certified): the pair with the same flip as the turn block sits at
Exp(−πU) with U/u = 2.12, 2.42, 2.68, 2.93, 2.99 at u = 0.01, 0.05, 0.1, 0.2, 0.3 and 3 at the strong end.
The pair's flip rises with the heat time (W2), so U is defined; U > 2u is W3, and W5 is U(2u) > 2·U(u).

## What it gives and what it does not

```text
the closed surface of the turn block      its diagonal step is exact: a boost by the pair's ratio (W1)          PROVED
(two dimensions, heat weight)             it is the pair's lost parts summed along the walk (W2)                PROVED
                                          it is the one law with seen L₀² and memory Σ2ⁿL_n² ; lost/seen < 1    PROVED
                                          flip positive, r < ρ < 1, at every coupling, from the strong side     PROVED
                                          the weak end is reached by walking back; no winding description (W4)  PROVED
                                          flip at least the mean step, at every coupling (W5 ; AG1-A8)          PROVED
                                          odd over even below its square under doubling (W6)                    PROVED
the fabric in four dimensions             NOT SHOWN. No exact block step is known there.
```

## Direction (a proposal, not a result)

AG1 asked for two inequalities per block step: source at most the arithmetic mean, flip at least the geometric
mean. For the one non-abelian record with exact data the second is now a theorem at every coupling, and the way
it was proved is the shape a proof in four dimensions would need: a certificate on the strong side (series of
one sign), an exact step, and a walk back whose loss is bounded once (W4's product) — no expansion at the weak
end. The compact part that blocks the expansion there (MG1's verdict; YM98 §5, YM99) is inside the walk. What
four dimensions lacks is the step: an inequality in place of W1 with a remainder of theorum/28's outward kind.
W1 says what to look for — the non-abelian record as the first change of an abelian one, stepping by a boost.

## What is not shown

- Nothing about the mass gap of the four-dimensional theory, and nothing about the Riemann hypothesis.
- W1 is the first change in the heat time of the pair's two step identities (Landen's), and W2 is, in older
  terms, Gauss's and Legendre's way to the second complete elliptic integral by the arithmetic–geometric mean
  (general knowledge; not re-read). The line's part: the turn block's closed surface as that object; the step
  as a boost by the pair's flip; the reading as lost parts along the diagonal walk; positivity and the weak end
  from the strong side without the mirror; and W5–W6, which I have not seen stated.
- W5 is a proof with exact arithmetic on a finite list of pieces (40 grid pieces, 100 boxes), not a term-by-term
  proof; the difference of its two sides is not of one sign in q.
- The limit π/2 in W4, the value at u = 1/16 and the 4u/π of W2 are agreement of enclosures or of floating
  numbers with the mirror description; equality is the classical seam law, pinned as in AG1.
- The torus record is a control at three couplings, not a theorem. U/u between 2 and 3 is floating; only U > 2u
  and U(2u) > 2U(u) are proved.
- The heat weight is used throughout. For the Wilson weight the faces do not compose exactly (CZ1-Z4).
- No mass count is derived (MG1).

## Claim boundary

```text
4a′K₊′ = aK₊ − bK₋ ; 4b′K₋′ = aK₋ − bK₊ ; ρ′r′ = (ρ − r)/(1 − rρ)             PROVED (lattice ; series to order 600)
16K₋ = b(L₀² + Σ2ⁿL_n²) ; 16K₊ = a(L₀² − Σ2ⁿL_n²) ; (K₋/b + K₊/a)/L² = 1/8     PROVED (from W1 ; series ; 60 digits at u = 1/4)
LOST/SEEN Λ < 1 ; 1 − Λ AGAINST 4u/π                                          PROVED at every coupling ; enclosures at four ; 4u/π is the mirror's
0 < r < ρ < 1 ; 3γ < α < γ + γ′ ; ρ(q) > F/S at q²                             PROVED at every coupling (walk from q ≤ ½)
β + 1 = 2(β′ + 1)/(1 + δ) ; WEAK END BY WALKING BACK                           PROVED ; agreement with the mirror to 80 digits
ρ′(1 + ρ) > 2√ρ                                                               PROVED at every coupling (exact pieces)
ρ(q²) > ρ(q) ; ε(q²) < ε(q)²                                                   PROVED at every coupling
TORUS RECORD: REVERSED                                                        CERTIFIED at q = 1/5, 1/3, 1/2 only
THE WALK IN FOUR DIMENSIONS ; A MASS GAP ; RH                                 NOT SHOWN ; OPEN
```

## Reproduce

```text
python tw1_turn_block_walk.py
python -m unittest test_tw1
```

## Later note (HL1, 9 October)

[HL1](../hl1/HL1_HELICAL_WALK.md) writes this stage in the terms of the UGD number and the EMK geometry. W1 is
the multiplicative law of the block K₊ + K₋K (EMK-2 T2; F00-G 7.1): Log x(ρ′r′) = Log x(ρ) − Log x(r). W2 is a
ledger over the sheets of the deck step 1 + ι, each with scale 2ⁿ and count L_n². The weak-end value 4u/π is
what a reading in powers of u keeps; its residue is certified there and is the first sheet term. The torus
record of the control is the two-turn commutator record. Nothing above is changed.

## Later note (SL1, 9 October)

[SL1](../sl1/SL1_SEAM_LAW_FROM_THE_WALK.md) proves the seam law of the pair — S(1/u) = u·S, F(1/u) = u·L, L(1/u) = u·F — from the diagonal
step, the two means and one count at the weak end; the only input from outside the walk is the area under
Exp(−πx²). So the limit π/2 of W4 and the 4u/π of W2 follow from the law and are no longer held as classical. Nothing above is changed.
